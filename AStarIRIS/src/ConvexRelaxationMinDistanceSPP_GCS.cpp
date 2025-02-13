#include "ConvexRelaxationMinDistanceSPP_GCS.h"
#include "Graph.h"
#include "PolyhedronNode.h"
#include "PolyhedronTerminalNode.h"
#include "PointNode.h"
#include "fusion.h"
#include <vector>
#include "EigenUtils.h"
#include <random>

ConvexRelaxationMinDistanceSPP_GCS::ConvexRelaxationMinDistanceSPP_GCS(Graph* g, const int &startKey, const int &targetKey): PerspectiveSPP_GCS(g,2,startKey,targetKey)
{
}

void ConvexRelaxationMinDistanceSPP_GCS::computeFeasibleSolution(const int& maxIters, const bool& simplifyGraph)
{
	int iters = 0;
	this->perspectiveSolution.cost = (double)std::numeric_limits<double>::infinity();
	Path_t OptimalPath;
	Path_t path;
	if (simplifyGraph)
		this->simplifyGraph();
	while (iters < maxIters)
	{
		if (!this->getMCPath(path))
		{
			//std::cout << "MCPath returned an invalid path" << std::endl;
			continue;
		}
		//int terminalNodeKey = path.nodeKeys[path.nodeKeys.size() - 2];
		//Node* node = this->g.getNode(terminalNodeKey);
		int nEdges = path.edgeKeys.size();
		int nVertex = path.nodeKeys.size();
		Model::t MPath = new Model("ConvexRelaxationMinDistance:FeasibleSolution");
		auto _M = finally([&]() { MPath->dispose(); });
		Variable::t lPath = MPath->variable("l", nEdges, Domain::greaterThan(0.0));
		Variable::t xPath = MPath->variable("x", Set::make(this->n, nVertex), Domain::unbounded());

		//std::cout << "Vertex" << std::endl;
		//std::cout << nVertex << std::endl;
		//std::cout << "Edges" << std::endl;
		//std::cout << nEdges << std::endl;
		//std::cout << "nodes keys ";
		//for (std::vector<int>::iterator it = this->optimalPath.nodeKeys.begin(); it != this->optimalPath.nodeKeys.end(); it++)
		//	std::cout << *it << " ";
		//std::cout << std::endl;
		//std::cout << "edge keys ";
		//for (std::vector<int>::iterator it = this->optimalPath.edgeKeys.begin(); it != this->optimalPath.edgeKeys.end(); it++)
		//	std::cout << *it << " ";
		//std::cout << std::endl;

		Expression::t x_expr = xPath->slice(new_array_ptr<int, 1>({ 0,0 }), new_array_ptr<int, 1>({ n,nVertex - 1 }))->transpose();
		Expression::t xp_expr = xPath->slice(new_array_ptr<int, 1>({ 0,1 }), new_array_ptr<int, 1>({ n,nVertex }))->transpose();
		MPath->constraint(Expr::hstack(lPath, Expr::sub(x_expr, xp_expr)), Domain::inQCone());

		for (int i = 0; i < nVertex; i++)
		{
			if ((path.nodeKeys[i] != startKey) && (path.nodeKeys[i] != targetKey))
			{
				PolyhedronNode* uNode = (PolyhedronNode*)this->g->getNode(path.nodeKeys[i]);
				Expression::t x_expr = xPath->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 }));
				std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(uNode->polyhedron.A);
				//std::shared_ptr<ndarray<double, 2>> b_ptr = Eigen2NdArray(Eigen::MatrixXd(uNode->polyhedron.b));
				std::shared_ptr<ndarray<double, 1>> b_ptr = Eigen2NdArray(uNode->polyhedron.b);
				MPath->constraint(Expr::sub(Expr::mul(A_ptr, x_expr), b_ptr), Domain::lessThan(0.0));
			}
			else if (path.nodeKeys[i] == startKey)
			{
				Expression::t x_expr = xPath->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 }));
				MPath->constraint(Expr::sub(x_expr, this->qstart), Domain::equalsTo(0.0));
			}
			else if (path.nodeKeys[i] == targetKey)
			{
				Expression::t x_expr = xPath->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 }));
				MPath->constraint(Expr::sub(x_expr, this->qtarget), Domain::equalsTo(0.0));
			}
		}
		MPath->objective(ObjectiveSense::Minimize, Expr::sum(lPath));
		//MPath->writeTask("dump.ptf");
		MPath->solve();
		this->feasibleSolution.status = MPath->getPrimalSolutionStatus();
		if (feasibleSolution.status == SolutionStatus::Optimal) {
			if (MPath->primalObjValue() < this->feasibleSolution.cost)
			{
				Eigen::VectorXd xvec = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(MPath->getVariable("x")->level()->begin(), this->n * nVertex));
				this->feasibleSolution.x = xvec.reshaped(nVertex, this->n).eval();
				this->feasibleSolution.cost = MPath->primalObjValue();
				this->optimalPath = path;
			}
		}
		else
		{
			std::cout << "Another solution status" << std::endl;
			std::cout << "Solution status: " << MPath->getPrimalSolutionStatus() << std::endl;
			//break;
		}		
		iters++;
	}
}

void ConvexRelaxationMinDistanceSPP_GCS::allocateVariables(Model::t& M)
{
	int nE = this->g->numEdges;
	this->l = M->variable("l", nE, Domain::greaterThan(0.0));
	this->z = M->variable("z", Set::make(this->n * this->N,nE), Domain::unbounded());
	this->y = M->variable("y", nE, Domain::greaterThan(0.0));
}

void ConvexRelaxationMinDistanceSPP_GCS::setConstraints(Model::t& M)
{
	int nE = this->g->numEdges;
	int nN = this->g->numNodes;
	M->constraint(Expr::hstack(l, Expr::sub(z->slice(new_array_ptr<int, 1>({ 0,0 }), new_array_ptr<int, 1>({ n,nE }))->transpose(), z->slice(new_array_ptr<int, 1>({ n,0 }), new_array_ptr<int, 1>({ 2 * n,nE }))->transpose())), Domain::inQCone());

	std::shared_ptr<ndarray<int, 1>>  idxsOut = new_array_ptr<int>(this->g->findOutEdges(this->startKey));
	std::shared_ptr<ndarray<int, 1>>  idxtIn = new_array_ptr<int>(this->g->findInEdges(this->targetKey));
	M->constraint(Expr::sum(y->pick(idxsOut)), Domain::equalsTo(1.0));
	M->constraint(Expr::sum(y->pick(idxtIn)), Domain::equalsTo(1.0));

	std::vector<int> nodeKeys = this->g->getNodeKeys();

	for (int i = 0; i < nN; i++)
	{
		if ((nodeKeys[i] != startKey) && (nodeKeys[i] != targetKey))
		{
			std::vector<int> outEdges = this->g->findOutEdges(nodeKeys[i]);
			std::vector<int> inEdges = this->g->findInEdges(nodeKeys[i]);
			std::shared_ptr<ndarray<int, 1>>  idxvOut = new_array_ptr<int>(outEdges);
			std::shared_ptr<ndarray<int, 1>>  idxvIn = new_array_ptr<int>(inEdges);
			M->constraint(Expr::sum(y->pick(idxvOut)), Domain::lessThan(1.0));
			Expression::t expr1(Expr::zeros(n));
			Expression::t expr2(Expr::zeros(n));
			for (int j = 0; j < inEdges.size(); j++)
			{
				expr1 = Expr::add(expr1, z->slice(new_array_ptr<int, 1>({ n,inEdges[j] }), new_array_ptr<int, 1>({ 2 * n,inEdges[j] + 1 })));
			}
			for (int j = 0; j < outEdges.size(); j++)
			{
				expr2 = Expr::add(expr2, z->slice(new_array_ptr<int, 1>({ 0,outEdges[j] }), new_array_ptr<int, 1>({ n,outEdges[j] + 1 })));
			}
			M->constraint(Expr::sub(Expr::sum(y->pick(idxvIn)), Expr::sum(y->pick(idxvOut))), Domain::equalsTo(0.0));
			M->constraint(Expr::sub(expr1, expr2), Domain::equalsTo(0.0));
		}
	}
	std::vector<Edge> edges(this->g->getEdges());
	std::vector<Edge>::iterator edgeIt = edges.begin();

	for (int i = 0; edgeIt != edges.end(); edgeIt++)
	{
		int u = edgeIt->first;
		int v = edgeIt->second;

		if (u == startKey)
		{
			Expression::t z_expr = z->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 }));
			Expression::t qy = Expr::mul(this->qstart, y->pick(new_array_ptr<int, 1>({ i })));
			M->constraint(Expr::sub(z_expr, qy), Domain::equalsTo(0.0));
		}
		else
		{
			PolyhedronNode* uNode = (PolyhedronNode*)this->g->getNode(u);
			Expression::t z_expr = z->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 }));
			std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(uNode->polyhedron.A);
			std::shared_ptr<ndarray<double, 2>> b_ptr = Eigen2NdArray(Eigen::MatrixXd(uNode->polyhedron.b));
			M->constraint(Expr::sub(Expr::mul(A_ptr, z_expr), Expr::mul(b_ptr, y->pick(new_array_ptr<int, 1>({ i })))), Domain::lessThan(0.0));
		}

		if (v == targetKey)
		{
			Expression::t p_expr = z->slice(new_array_ptr<int, 1>({ n,i }), new_array_ptr<int, 1>({ 2 * n,i + 1 }));
			Expression::t qy = Expr::mul(this->qtarget, y->pick(new_array_ptr<int, 1>({ i })));
			M->constraint(Expr::sub(p_expr, qy), Domain::equalsTo(0.0));
		}
		else
		{
			PolyhedronNode* vNode = (PolyhedronNode*)this->g->getNode(v);
			Expression::t p_expr = z->slice(new_array_ptr<int, 1>({ n,i }), new_array_ptr<int, 1>({ 2 * n,i + 1 }));
			std::shared_ptr<ndarray<double, 2>> Ap_ptr = Eigen2NdArray(vNode->polyhedron.A);
			std::shared_ptr<ndarray<double, 2>> bp_ptr = Eigen2NdArray(Eigen::MatrixXd(vNode->polyhedron.b));
			M->constraint(Expr::sub(Expr::mul(Ap_ptr, p_expr), Expr::mul(bp_ptr, y->pick(new_array_ptr<int, 1>({ i })))), Domain::lessThan(0.0));
		}
		i++;
	}
}

void ConvexRelaxationMinDistanceSPP_GCS::setObjective(Model::t& M)
{
	M->objective(ObjectiveSense::Minimize, Expr::add(Expr::sum(this->l), Expr::mul(Expr::sum(this->y), 1e-4)));
}

void ConvexRelaxationMinDistanceSPP_GCS::optimize()
{
		Model::t M = new Model("ConvexRelaxationMinDistance");
		auto _M = finally([&]() { M->dispose(); });
		int nE = this->g->numEdges;
		
		this->allocateVariables(M);
		this->setConstraints(M);
		this->setObjective(M);
		M->solve();
		this->perspectiveSolution.status = M->getPrimalSolutionStatus();
		if (perspectiveSolution.status == SolutionStatus::Optimal) {
			Eigen::VectorXd y = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(M->getVariable("y")->level()->begin(), nE));
			Eigen::VectorXd l = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(M->getVariable("l")->level()->begin(), nE));
			Eigen::VectorXd zvec = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(M->getVariable("z")->level()->begin(), this->n * nE * this->N));
			this->perspectiveSolution.y = y;
			this->perspectiveSolution.l = l;
			this->perspectiveSolution.z = zvec.reshaped(nE, this->n * this->N).eval();
			this->perspectiveSolution.N = this->N;
			this->perspectiveSolution.cost = M->primalObjValue();
			this->perspectiveSolutionSolved = true;
		}
		else
		{
			std::cout << "Another solution status" << std::endl;
			std::cout << "Solution status: " << M->getPrimalSolutionStatus() << std::endl;
		}
		M->dispose();
}
