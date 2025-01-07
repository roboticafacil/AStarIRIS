#include "MinDistanceSolver.h"
#include "Graph.h"
#include "PolyhedronNode.h"
#include "PointNode.h"
#include "fusion.h"
#include <vector>
#include "EigenNdArray.h"

MinDistanceSolver::MinDistanceSolver(const Graph& g, const int& startKey, const int& targetKey) : Solver(g, startKey, targetKey)
{
}

void MinDistanceSolver::setTask()
{
	this->M = new Model("MinDistance");
	this->solverAllocated = true;
	this->l = this->M->variable("l", this->nE, Domain::greaterThan(0.0));
	this->z = this->M->variable("z", Set::make(this->n, this->nE), Domain::unbounded());
	this->p = this->M->variable("p", Set::make(this->n, this->nE), Domain::unbounded());
	this->y = this->M->variable("y", this->nE, Domain::integral(Domain::inRange(0,1)));
	//this->y = this->M->variable("y", this->nE, Domain::binary());

	this->M->constraint(Expr::hstack(l, Expr::sub(z->transpose(), p->transpose())), Domain::inQCone());


	std::shared_ptr<ndarray<int, 1>>  idxsOut = new_array_ptr<int>(this->g.findOutEdges(this->startKey));
	std::shared_ptr<ndarray<int, 1>>  idxtIn = new_array_ptr<int>(this->g.findInEdges(this->targetKey));
	this->M->constraint(Expr::sum(y->pick(idxsOut)), Domain::equalsTo(1));
	this->M->constraint(Expr::sum(y->pick(idxtIn)), Domain::equalsTo(1));

	std::vector<int> nodeKeys = this->g.getNodeKeys();

	for (int i = 0; i < this->nN; i++)
	{
		if ((nodeKeys[i] != startKey) && (nodeKeys[i] != targetKey))
		{
			std::vector<int> outEdges = this->g.findOutEdges(nodeKeys[i]);
			std::vector<int> inEdges = this->g.findInEdges(nodeKeys[i]);
			std::shared_ptr<ndarray<int, 1>>  idxvOut = new_array_ptr<int>(outEdges);
			std::shared_ptr<ndarray<int, 1>>  idxvIn = new_array_ptr<int>(inEdges);
			this->M->constraint(Expr::sum(y->pick(idxvOut)), Domain::lessThan(1.0));
			Expression::t expr1(Expr::zeros(n));
			Expression::t expr2(Expr::zeros(n));
			for (int j = 0; j < inEdges.size(); j++)
			{
				expr1 = Expr::add(expr1, p->slice(new_array_ptr<int, 1>({ 0,inEdges[j] }), new_array_ptr<int, 1>({ n,inEdges[j] + 1 })));
			}
			for (int j = 0; j < outEdges.size(); j++)
			{
				expr2 = Expr::add(expr2, z->slice(new_array_ptr<int, 1>({ 0,outEdges[j] }), new_array_ptr<int, 1>({ n,outEdges[j] + 1 })));
			}
			this->M->constraint(Expr::sub(Expr::sum(y->pick(idxvIn)), Expr::sum(y->pick(idxvOut))), Domain::equalsTo(0));
			this->M->constraint(Expr::sub(expr1, expr2), Domain::equalsTo(0.0));
		}
	}
	std::vector<Edge> edges(this->g.getEdges());
	std::vector<Edge>::iterator edgeIt = edges.begin();

	for (int i = 0; edgeIt != edges.end(); edgeIt++)
	{
		int u = edgeIt->first;
		int v = edgeIt->second;

		if (u == startKey)
		{
			Expression::t z_expr = z->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 }));
			Expression::t qy = Expr::mul(this->qstart, y->pick(new_array_ptr<int, 1>({ i })));
			this->M->constraint(Expr::sub(z_expr, qy), Domain::equalsTo(0.0));
		}
		else
		{
			PolyhedronNode* uNode = (PolyhedronNode*)this->g.getNode(u);
			Expression::t z_expr = z->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 }));
			std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(uNode->polyhedron.A);
			std::shared_ptr<ndarray<double, 2>> b_ptr = Eigen2NdArray(Eigen::MatrixXd(uNode->polyhedron.b));
			this->M->constraint(Expr::sub(Expr::mul(A_ptr, z_expr), Expr::mul(b_ptr, y->pick(new_array_ptr<int, 1>({ i })))), Domain::lessThan(0.0));
		}

		if (v == targetKey)
		{
			Expression::t p_expr = p->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 }));
			Expression::t qy = Expr::mul(this->qtarget, y->pick(new_array_ptr<int, 1>({ i })));
			this->M->constraint(Expr::sub(p_expr, qy), Domain::equalsTo(0.0));
		}
		else
		{
			PolyhedronNode* vNode = (PolyhedronNode*)this->g.getNode(v);
			Expression::t p_expr = p->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 }));
			std::shared_ptr<ndarray<double, 2>> Ap_ptr = Eigen2NdArray(vNode->polyhedron.A);
			std::shared_ptr<ndarray<double, 2>> bp_ptr = Eigen2NdArray(Eigen::MatrixXd(vNode->polyhedron.b));
			this->M->constraint(Expr::sub(Expr::mul(Ap_ptr, p_expr), Expr::mul(bp_ptr, y->pick(new_array_ptr<int, 1>({ i })))), Domain::lessThan(0.0));
		}
		i++;
	}
	this->M->objective(ObjectiveSense::Minimize, Expr::sum(this->l));
	//this->M->writeTask("dump.ptf");
}

void MinDistanceSolver::computeFeasibleSolution(const int &maxIters)
{
	int nextKey = startKey;

	while (nextKey != targetKey)
	{
		this->optimalPath.nodeKeys.push_back(nextKey);
		std::vector<int> outEdges = this->g.findOutEdges(nextKey);
		double max_y = 0.;
		int best_edge = -1;
		for (std::vector<int>::iterator it = outEdges.begin(); it != outEdges.end(); it++)
		{
			double current_y = this->MIPSolution.y(*it);
			if (current_y > max_y)
			{
				max_y = current_y;
				best_edge = *it;
			}
		}
		nextKey = this->g.getEdge(best_edge).second;
		this->optimalPath.edgeKeys.push_back(best_edge);
	}
	this->optimalPath.nodeKeys.push_back(nextKey);

	this->feasibleSolution.x.resize(this->optimalPath.nodeKeys.size(), this->MIPSolution.x.cols());
	int i = 0;
	this->feasibleSolution.cost = 0.;
	std::vector<int> nodeKeys = this->g.getNodeKeys();
	for (std::vector<int>::iterator it = this->optimalPath.nodeKeys.begin(); it != this->optimalPath.nodeKeys.end(); it++)
	{
		std::vector<int>::iterator it1 = std::find(nodeKeys.begin(), nodeKeys.end(), *it);
		int j = it1- nodeKeys.begin();
		this->feasibleSolution.x.row(i) = this->MIPSolution.x.row(j);
		i++;
	}
	i = 0;
	for (std::vector<int>::iterator it = this->optimalPath.edgeKeys.begin(); it != this->optimalPath.edgeKeys.end(); it++)
	{
		int j = *it;
		this->feasibleSolution.cost += this->MIPSolution.l(j);
		i++;
	}
}

void MinDistanceSolver::solve()
{
	this->M->solve();
	this->MIPSolution.status = this->M->getPrimalSolutionStatus();
	if (this->MIPSolution.status == SolutionStatus::Optimal) {
		Eigen::VectorXd y = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>((this->M->getVariable("y")->level()->begin()), this->nE));
		Eigen::VectorXd l = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(this->M->getVariable("l")->level()->begin(), this->nE));
		Eigen::VectorXd zvec = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(this->M->getVariable("z")->level()->begin(), this->n * this->nE));
		Eigen::VectorXd pvec = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(this->M->getVariable("p")->level()->begin(), this->n * this->nE));
		this->MIPSolution.y = y;
		this->MIPSolution.l = l;
		this->MIPSolution.z = zvec.reshaped(this->nE, this->n).eval();
		this->MIPSolution.p = pvec.reshaped(this->nE, this->n).eval();
		this->MIPSolution.x = Eigen::MatrixXd::Zero(this->nN, this->n);
		std::vector<int> nodeKeys = this->g.getNodeKeys();
		for (int i = 0; i < this->nN; i++)
		{
			std::vector<int> inEdges = this->g.findInEdges(nodeKeys[i]);
			if (inEdges.size() > 0)
			{
				for (std::vector<int>::iterator it = inEdges.begin(); it != inEdges.end(); it++)
				{
					this->MIPSolution.x.row(i) += this->MIPSolution.y(*it) * this->MIPSolution.p.row(*it);
				}
			}
			else
			{
				std::vector<int> outEdges = this->g.findOutEdges(nodeKeys[i]);
				if (outEdges.size() > 0)
				{
					for (std::vector<int>::iterator it = outEdges.begin(); it != outEdges.end(); it++)
					{
						this->MIPSolution.x.row(i) += this->MIPSolution.y(*it) * this->MIPSolution.z.row(*it);
					}
				}
				else
				{
					std::cout << "Node " << i << " has no edges" << std::endl;
				}
			}
		}
		this->MIPSolution.cost = this->M->primalObjValue();
	}
	else
	{
		std::cout << "Another solution status" << std::endl;
		std::cout << "Solution status: " << this->M->getPrimalSolutionStatus() << std::endl;
	}
	this->computeFeasibleSolution();
}