#include <Eigen/Cholesky>
#include "ConvexRelaxationBezierCurveSPP_GCS.h"
#include "Graph.h"
#include "PolyhedronNode.h"
#include "PolyhedronTerminalNode.h"
#include "PointNode.h"
#include "Range.h"
#include "fusion.h"
#include <vector>
#include "EigenUtils.h"
#include <random>

ConvexRelaxationBezierCurveSPP_GCS::ConvexRelaxationBezierCurveSPP_GCS(Graph* navGraph, Graph* gcs, std::map<std::pair<int, int>, int>* navGraph2gcs, const ConvexRelaxationBezierCurveParams_t& params, const int& startKey, const int& targetKey, const Eigen::MatrixXd& CpStart, const Eigen::MatrixXd& CpTarget) : PerspectiveSPP_GCS(navGraph,params.N,startKey, targetKey), gcs(gcs), navGraph2gcs(navGraph2gcs), params(params)
{
	if (CpStart.rows() >= 1)
	{
		std::cout << Eigen::VectorXd(CpStart.row(0)) << std::endl;
		cpStart = Eigen2NdArray(Eigen::VectorXd(CpStart.row(0)));
	}
	else
		cpStart = Eigen2NdArray(Eigen::VectorXd(0));
	if (CpStart.rows() >= 1)
		cpTarget=Eigen2NdArray(Eigen::VectorXd(CpTarget.row(0)));
	else
		cpTarget = Eigen2NdArray(Eigen::VectorXd(0));
	if (CpStart.rows() >= 2)
		cppStart = Eigen2NdArray(Eigen::VectorXd(CpStart.row(1)));
	else
		cppStart = Eigen2NdArray(Eigen::VectorXd(0));
	if (CpStart.rows() >= 2)
		cppTarget = Eigen2NdArray(Eigen::VectorXd(CpTarget.row(1)));
	else
		cppTarget = Eigen2NdArray(Eigen::VectorXd(0));
	if (CpStart.rows() >= 3)
		cpppStart = Eigen2NdArray(Eigen::VectorXd(CpStart.row(2)));
	else
		cpppStart = Eigen2NdArray(Eigen::VectorXd(0));
	if (CpStart.rows() >= 3)
		cpppTarget = Eigen2NdArray(Eigen::VectorXd(CpTarget.row(2)));
	else
		cpppTarget = Eigen2NdArray(Eigen::VectorXd(0));
	this->computeBezierL2NormMatrix();
}

void ConvexRelaxationBezierCurveSPP_GCS::computeFeasibleSolution(const int& maxIters, const bool& simplifyGraph)
{
	int iters = 0;
	this->feasibleSolution.cost = (double)std::numeric_limits<double>::infinity();
	Path_t OptimalPath;
	Path_t path;
	while (iters < maxIters)
	{
		if (!this->getMCPath(path))
		{
			//std::cout << "MCPath returned an invalid path" << std::endl;
			continue;
		}
		int nEdges = path.edgeKeys.size();
		int nVertex = path.nodeKeys.size();
		Model::t MPath = new Model("ConvexRelaxationBezierCurve:FeasibleSolution");
		auto _M = finally([&]() { MPath->dispose(); });
		Variable::t lPath = MPath->variable("l", nEdges, Domain::greaterThan(0.0));
		Variable::t xPath = MPath->variable("x", Set::make(this->n*this->N, nEdges), Domain::unbounded());

		assert((this->N - this->params.costDerivativeOrder) >= 2);
		assert((this->N - this->params.continuityDerivativeOrder) >= 2);

		if (this->params.costDerivativeOrder == 1)
		{
			Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({ n,0 }), new_array_ptr<int, 1>({ this->N * n,nEdges }))->transpose();
			Variable::t x0 = xPath->slice(new_array_ptr<int, 1>({ 0,0 }), new_array_ptr<int, 1>({ (this->N - 1) * n,nEdges }))->transpose();
			//this->y->pick()
			
			if (this->params.squaredRootL2NormCost)
				MPath->constraint(Expr::hstack(lPath, Expr::mul(Expr::sub(x1, x0), Qptr)), Domain::inQCone());
			else
			{
				//MPath->constraint(Expr::hstack(Expr::constTerm(nEdges, 1.), lPath, Expr::mul(Expr::sub(x1, x0), Qptr)), Domain::inRotatedQCone());
				MPath->constraint(Expr::hstack(y, lPath, Expr::mul(Expr::sub(x1, x0), Qptr)), Domain::inRotatedQCone());
			}
		}
		else if (this->params.costDerivativeOrder == 2)
		{
			Variable::t x2 = xPath->slice(new_array_ptr<int, 1>({ 2 * n,0 }), new_array_ptr<int, 1>({ this->N * n,nEdges }))->transpose();
			Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({ n,0 }), new_array_ptr<int, 1>({ (this->N - 1) * n,nEdges }))->transpose();
			Variable::t x0 = xPath->slice(new_array_ptr<int, 1>({ 0,0 }), new_array_ptr<int, 1>({ (this->N - 2) * n,nEdges }))->transpose();
			if (this->params.squaredRootL2NormCost)
				MPath->constraint(Expr::hstack(lPath, Expr::mul(Expr::add(Expr::sub(x2, Expr::mul(2, x1)), x0), Qptr)), Domain::inQCone());
			else
			{
				//MPath->constraint(Expr::hstack(Expr::constTerm(nEdges, 1.), lPath, Expr::mul(Expr::add(Expr::sub(x2, Expr::mul(2, x1)), x0), Qptr)), Domain::inRotatedQCone());
				MPath->constraint(Expr::hstack(y, lPath, Expr::mul(Expr::add(Expr::sub(x2, Expr::mul(2, x1)), x0), Qptr)), Domain::inRotatedQCone());
			}
		}
		else if (this->params.costDerivativeOrder == 3)
		{
			Variable::t x3 = xPath->slice(new_array_ptr<int, 1>({ 3 * n,0 }), new_array_ptr<int, 1>({ this->N * n,nEdges }))->transpose();
			Variable::t x2 = xPath->slice(new_array_ptr<int, 1>({ 2 * n,0 }), new_array_ptr<int, 1>({ (this->N - 1) * n,nEdges }))->transpose();
			Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({ n,0 }), new_array_ptr<int, 1>({ (this->N - 2) * n,nEdges }))->transpose();
			Variable::t x0 = xPath->slice(new_array_ptr<int, 1>({ 0,0 }), new_array_ptr<int, 1>({ (this->N - 3) * n,nEdges }))->transpose();
			if (this->params.squaredRootL2NormCost)
				MPath->constraint(Expr::hstack(lPath, Expr::mul(Expr::add(Expr::sub(x3, Expr::mul(3, x2)), Expr::sub(Expr::mul(3, x1), x0)), Qptr)), Domain::inQCone());
			else
			{
				//MPath->constraint(Expr::hstack(Expr::constTerm(nEdges, 1.), lPath, Expr::mul(Expr::add(Expr::sub(x3, Expr::mul(3, x2)), Expr::sub(Expr::mul(3, x1), x0)), Qptr)), Domain::inRotatedQCone());
				MPath->constraint(Expr::hstack(y, lPath, Expr::mul(Expr::add(Expr::sub(x3, Expr::mul(3, x2)), Expr::sub(Expr::mul(3, x1), x0)), Qptr)), Domain::inRotatedQCone());
			}
		}

		Variable::t apS, apT, appS, appT, apppS, apppT;
		if ((cpStart->size() > 0) && (this->params.continuityDerivativeOrder >= 1))
		{
			MPath->variable("apS", 1, Domain::greaterThan(0.0));
		}
		if ((cpTarget->size() > 0) && (this->params.continuityDerivativeOrder >= 1))
		{
			MPath->variable("apT", 1, Domain::greaterThan(0.0));
		}
		if ((cppStart->size() > 0) && (this->params.continuityDerivativeOrder >= 2))
		{
			MPath->variable("appS", 1, Domain::greaterThan(0.0));
		}
		if ((cppTarget->size() > 0) && (this->params.continuityDerivativeOrder >= 2))
		{
			MPath->variable("appT", 1, Domain::greaterThan(0.0));
		}
		if ((cpppStart->size() > 0) && (this->params.continuityDerivativeOrder >= 3))
		{
			MPath->variable("apppS", 1, Domain::greaterThan(0.0));
		}
		if ((cpppTarget->size() > 0) && (this->params.continuityDerivativeOrder >= 3))
		{
			MPath->variable("apppT", 1, Domain::greaterThan(0.0));
		}

		//StartNode
		if ((cpStart->size() > 0)&& (this->params.continuityDerivativeOrder >= 1))
		{
			Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({ n,0 }), new_array_ptr<int, 1>({ 2 * n,1 }));
			Variable::t x0 = xPath->slice(new_array_ptr<int, 1>({ 0,0 }), new_array_ptr<int, 1>({ n,1 }));
			MPath->constraint(Expr::sub(Expr::sub(x1, x0), Expr::mul(apS, cpStart)), Domain::equalsTo(0.0));
		}
		if ((cppStart->size() > 0) && (this->params.continuityDerivativeOrder >= 2))
		{
			Variable::t x2 = xPath->slice(new_array_ptr<int, 1>({ 2 * n,0}), new_array_ptr<int, 1>({ 3 * n,1 }));
			Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({ n,0}), new_array_ptr<int, 1>({ 2 * n,1 }));
			Variable::t x0 = xPath->slice(new_array_ptr<int, 1>({ 0,0}), new_array_ptr<int, 1>({ n,1 }));
			MPath->constraint(Expr::sub(Expr::add(Expr::sub(x2, Expr::mul(2, x1)), x0), Expr::mul(appS, cppStart)), Domain::equalsTo(0.0));
		}
		if ((cpppStart->size() > 0) && (this->params.continuityDerivativeOrder >= 3))
		{
			Variable::t x3 = xPath->slice(new_array_ptr<int, 1>({ 3 * n,0}), new_array_ptr<int, 1>({ 4 * n,1 }));
			Variable::t x2 = xPath->slice(new_array_ptr<int, 1>({ 2 * n,0}), new_array_ptr<int, 1>({ 3 * n,1 }));
			Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({ n,0}), new_array_ptr<int, 1>({ 2 * n,1 }));
			Variable::t x0 = xPath->slice(new_array_ptr<int, 1>({ 0,0}), new_array_ptr<int, 1>({ n,1 }));
			MPath->constraint(Expr::sub(Expr::add(Expr::sub(x3, Expr::mul(3, x2)), Expr::sub(Expr::mul(3, x1), x0)), Expr::mul(apppS, cpppStart)), Domain::equalsTo(0.0));
		}
		//Nodes in the middle of the path
		for (int i = 1; i <= (nEdges-1); i++)
		{
			//Last control point of previous edge
			Variable::t x0 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 1) * n,i-1}), new_array_ptr<int, 1>({ this->N * n,i}));
			//First control point of this edge
			Variable::t x0p = xPath->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 })); 
			MPath->constraint(Expr::sub(x0, x0p), Domain::equalsTo(0.0));
			if (this->params.continuityDerivativeOrder >= 1)
			{
				//First-order derivatives for last control points of previous edge
				Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 2) * n,i-1}), new_array_ptr<int, 1>({ (this->N - 1) * n,i}));
				//First-order derivatives for first control points of this edge
				Variable::t x1p = xPath->slice(new_array_ptr<int, 1>({ n,i }), new_array_ptr<int, 1>({ 2 * n,i + 1 }));
				MPath->constraint(Expr::sub(Expr::sub(x0, x1), Expr::sub(x1p, x0p)), Domain::equalsTo(0.0));
			}
			if (this->params.continuityDerivativeOrder >= 2)
			{
				//Second-order derivatives for last control points
				Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 2) * n,i-1}), new_array_ptr<int, 1>({ (this->N - 1) * n,i}));
				Variable::t x2 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 3) * n,i-1}), new_array_ptr<int, 1>({ (this->N - 2) * n,i}));
				Variable::t x1p = xPath->slice(new_array_ptr<int, 1>({ n,i }), new_array_ptr<int, 1>({ 2 * n,i + 1 }));
				Variable::t x2p = xPath->slice(new_array_ptr<int, 1>({ 2 * n,i }), new_array_ptr<int, 1>({ 3 * n,i + 1 }));
				MPath->constraint(Expr::sub(Expr::add(Expr::sub(x2, Expr::mul(2, x1)), x0), Expr::add(Expr::sub(x2p, Expr::mul(2, x1p)), x0p)), Domain::equalsTo(0.0));
			}
			if (this->params.continuityDerivativeOrder >= 3)
			{
				//Third-order derivatives for last control points
				Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 2) * n,i-1}), new_array_ptr<int, 1>({ (this->N - 1) * n,i}));
				Variable::t x2 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 3) * n,i-1}), new_array_ptr<int, 1>({ (this->N - 2) * n,i}));
				Variable::t x3 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 4) * n,i-1}), new_array_ptr<int, 1>({ (this->N - 4) * n,i}));
				Variable::t x1p = xPath->slice(new_array_ptr<int, 1>({ n,i }), new_array_ptr<int, 1>({ 2 * n,i + 1 }));
				Variable::t x2p = xPath->slice(new_array_ptr<int, 1>({ 2 * n,i }), new_array_ptr<int, 1>({ 3 * n,i + 1 }));
				Variable::t x3p = xPath->slice(new_array_ptr<int, 1>({ 3 * n,i }), new_array_ptr<int, 1>({ 4 * n,i + 1 }));
				MPath->constraint(Expr::sub(Expr::add(Expr::sub(x3, Expr::mul(3, x2)), Expr::sub(Expr::mul(3, x1), x0)), Expr::add(Expr::sub(x3p, Expr::mul(3, x2p)), Expr::sub(Expr::mul(3, x1p), x0p))), Domain::equalsTo(0.0));
			}
		}

		//Target node
		if ((cpTarget->size() > 0) && (this->params.continuityDerivativeOrder >= 1))
		{
			Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({(this->N-1)*n,nEdges-1}), new_array_ptr<int, 1>({ this->N * n,nEdges}));
			Variable::t x0 = xPath->slice(new_array_ptr<int, 1>({(this->N-2)*n,nEdges - 1 }), new_array_ptr<int, 1>({ (this->N - 1)*n,nEdges}));
			MPath->constraint(Expr::sub(Expr::sub(x1, x0), Expr::mul(apT, cpTarget)), Domain::equalsTo(0.0));
		}
		if ((cppTarget->size() > 0) && (this->params.continuityDerivativeOrder >= 2))
		{
			Variable::t x2 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 1) * n,nEdges - 1 }), new_array_ptr<int, 1>({ this->N * n,nEdges }));
			Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 2) * n,nEdges - 1 }), new_array_ptr<int, 1>({ (this->N - 1) * n * n,nEdges }));
			Variable::t x0 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 3) * n,nEdges - 1 }), new_array_ptr<int, 1>({ (this->N - 2) * n,nEdges }));
			MPath->constraint(Expr::sub(Expr::add(Expr::sub(x2, Expr::mul(2, x1)), x0), Expr::mul(appT, cppTarget)), Domain::equalsTo(0.0));
		}
		if ((cpppTarget->size() > 0) && (this->params.continuityDerivativeOrder >= 3))
		{
			Variable::t x3 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 1) * n,nEdges - 1 }), new_array_ptr<int, 1>({ this->N * n,nEdges }));
			Variable::t x2 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 2) * n,nEdges - 1 }), new_array_ptr<int, 1>({ (this->N - 1)* n,nEdges }));
			Variable::t x1 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 3) * n,nEdges - 1 }), new_array_ptr<int, 1>({ (this->N - 2) * n,nEdges }));
			Variable::t x0 = xPath->slice(new_array_ptr<int, 1>({ (this->N - 4) * n,nEdges - 1 }), new_array_ptr<int, 1>({ (this->N - 3) * n,nEdges }));
			MPath->constraint(Expr::sub(Expr::add(Expr::sub(x3, Expr::mul(3, x2)), Expr::sub(Expr::mul(3, x1), x0)), Expr::mul(apppT, cpppTarget)), Domain::equalsTo(0.0));
		}

		int i = 0;
		for (std::vector<int>::iterator edgeIt = path.edgeKeys.begin(); edgeIt != path.edgeKeys.end(); edgeIt++)
		{
			Edge edge=this->g->getEdge(*edgeIt);
			int u = edge.first;
			int v = edge.second;

			//First control point
			if (u == startKey)
			{
				Expression::t x_expr = xPath->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 }));
				MPath->constraint(Expr::sub(x_expr, this->qstart), Domain::equalsTo(0.0));
			}
			else
			{
				PolyhedronNode* uNode = (PolyhedronNode*)this->g->getNode(u);
				Expression::t x_expr = xPath->slice(new_array_ptr<int, 1>({ 0,i }), new_array_ptr<int, 1>({ n,i + 1 }));
				std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(uNode->polyhedron.A);
				std::shared_ptr<ndarray<double, 2>> b_ptr = Eigen2NdArray(Eigen::MatrixXd(uNode->polyhedron.b));
				MPath->constraint(Expr::sub(Expr::mul(A_ptr, x_expr), b_ptr), Domain::lessThan(0.0));
			}

			//Control points in the middle
			std::map<std::pair<int, int>, int>::iterator itNavGraph2GCS = this->navGraph2gcs->find(edge);
			if (itNavGraph2GCS != this->navGraph2gcs->end())
			{
				int gcsKey = itNavGraph2GCS->second;
				PolyhedronNode* node = (PolyhedronNode*)this->gcs->getNode(gcsKey);
				std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(node->polyhedron.A);
				std::shared_ptr<ndarray<double, 2>> b_ptr = Eigen2NdArray(Eigen::MatrixXd(node->polyhedron.b));
				for (int j = 1; j < (N - 1); j++)
				{
					Expression::t x_expr = xPath->slice(new_array_ptr<int, 1>({ j * n,i }), new_array_ptr<int, 1>({ (j + 1) * n,i + 1 }));
					MPath->constraint(Expr::sub(Expr::mul(A_ptr, x_expr), b_ptr), Domain::lessThan(0.0));
				}
			}
			else
			{
				Range* range = Range::getInstance();
				Eigen::MatrixXd ARange;
				Eigen::VectorXd bRange;
				range->getConstraints(ARange, bRange);
				std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(ARange);
				std::shared_ptr<ndarray<double, 2>> b_ptr = Eigen2NdArray(Eigen::MatrixXd(bRange));
				for (int j = 1; j < (N - 1); j++)
				{
					Expression::t x_expr = xPath->slice(new_array_ptr<int, 1>({ j * n,i }), new_array_ptr<int, 1>({ (j + 1) * n,i + 1 }));
					MPath->constraint(Expr::sub(Expr::mul(A_ptr, x_expr), b_ptr), Domain::lessThan(0.0));
				}
			}

			//Last control point
			if (v == targetKey)
			{
				Expression::t x_expr = xPath->slice(new_array_ptr<int, 1>({ (this->N - 1) * n,i }), new_array_ptr<int, 1>({ this->N * n,i + 1 }));
				MPath->constraint(Expr::sub(x_expr,this->qtarget), Domain::equalsTo(0.0));
			}
			else
			{
				PolyhedronNode* vNode = (PolyhedronNode*)this->g->getNode(v);
				Expression::t x_expr = xPath->slice(new_array_ptr<int, 1>({ (this->N - 1) * n,i}), new_array_ptr<int, 1>({ this->N * n,i + 1 }));
				std::shared_ptr<ndarray<double, 2>> Ap_ptr = Eigen2NdArray(vNode->polyhedron.A);
				std::shared_ptr<ndarray<double, 2>> bp_ptr = Eigen2NdArray(Eigen::MatrixXd(vNode->polyhedron.b));
				MPath->constraint(Expr::sub(Expr::mul(Ap_ptr, x_expr), bp_ptr), Domain::lessThan(0.0));
			}
			i++;
		}


		MPath->objective(ObjectiveSense::Minimize, Expr::sum(lPath));
		//MPath->writeTask("dump.ptf");
		MPath->solve();
		this->feasibleSolution.status = MPath->getPrimalSolutionStatus();
		if (feasibleSolution.status == SolutionStatus::Optimal) {
			if (MPath->primalObjValue() < this->feasibleSolution.cost)
			{
				//Eigen::VectorXd xvec = Eigen::VectorXd(Eigen::Map<Eigen::VectorXd>(MPath->getVariable("x")->level()->begin(), this->n * this->N * nEdges));
				Eigen::MatrixXd xvec = Eigen::MatrixXd(Eigen::Map<Eigen::MatrixXd>(MPath->getVariable("x")->level()->begin(), nEdges, this->n * this->N));
				xvec=xvec.transpose().eval().reshaped(n,this->N*nEdges).eval().transpose().eval();
				this->feasibleSolution.x = xvec;
				//std::cout << this->feasibleSolution.x << std::endl;
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

void ConvexRelaxationBezierCurveSPP_GCS::allocateVariables(Model::t& M)
{
	int nE = this->g->numEdges;
	this->l = M->variable("l", nE, Domain::greaterThan(0.0));
	this->z = M->variable("z", Set::make(this->n * this->N, nE), Domain::unbounded());
	this->y = M->variable("y", nE, Domain::greaterThan(0.0));
}

void ConvexRelaxationBezierCurveSPP_GCS::setConstraints(Model::t& M)
{
	int nE = this->g->numEdges;
	int nN = this->g->numNodes;
	assert((this->N - this->params.costDerivativeOrder) >= 2);
	assert((this->N - this->params.continuityDerivativeOrder) >= 2);


	if (this->params.costDerivativeOrder == 1)
	{
		Variable::t z1 = this->z->slice(new_array_ptr<int, 1>({ n,0 }), new_array_ptr<int, 1>({ this->N * n,nE }))->transpose();
		Variable::t z0 = this->z->slice(new_array_ptr<int, 1>({ 0,0 }), new_array_ptr<int, 1>({ (this->N - 1) * n,nE }))->transpose();
		if (this->params.squaredRootL2NormCost)
		{
			M->constraint(Expr::hstack(this->l, Expr::mul(Expr::sub(z1, z0), Qptr)), Domain::inQCone());
		}
		else
		{
			M->constraint(Expr::hstack(Expr::constTerm(nE, 1.), this->l, Expr::mul(Expr::sub(z1, z0), Qptr)), Domain::inRotatedQCone());
		}
	}
	else if (this->params.costDerivativeOrder == 2)
	{
		Variable::t z2 = this->z->slice(new_array_ptr<int, 1>({ 2 * n,0 }), new_array_ptr<int, 1>({ this->N * n,nE }))->transpose();
		Variable::t z1 = this->z->slice(new_array_ptr<int, 1>({ n,0 }), new_array_ptr<int, 1>({ (this->N - 1) * n,nE }))->transpose();
		Variable::t z0 = this->z->slice(new_array_ptr<int, 1>({ 0,0 }), new_array_ptr<int, 1>({ (this->N - 2) * n,nE }))->transpose();
		if (this->params.squaredRootL2NormCost)
		{
			M->constraint(Expr::hstack(this->l, Expr::mul(Expr::add(Expr::sub(z2, Expr::mul(2, z1)), z0), Qptr)), Domain::inQCone());
		}
		else
		{
			M->constraint(Expr::hstack(Expr::constTerm(nE, 1.), this->l, Expr::mul(Expr::add(Expr::sub(z2, Expr::mul(2, z1)), z0), Qptr)), Domain::inRotatedQCone());
		}
	}
	else if (this->params.costDerivativeOrder == 3)
	{
		Variable::t z3 = this->z->slice(new_array_ptr<int, 1>({ 3 * n,0 }), new_array_ptr<int, 1>({ this->N * n,nE }))->transpose();
		Variable::t z2 = this->z->slice(new_array_ptr<int, 1>({ 2 * n,0 }), new_array_ptr<int, 1>({ (this->N - 1) * n,nE }))->transpose();
		Variable::t z1 = this->z->slice(new_array_ptr<int, 1>({ n,0 }), new_array_ptr<int, 1>({ (this->N - 2) * n,nE }))->transpose();
		Variable::t z0 = this->z->slice(new_array_ptr<int, 1>({ 0,0 }), new_array_ptr<int, 1>({ (this->N - 3) * n,nE }))->transpose();
		if (this->params.squaredRootL2NormCost)
		{
			M->constraint(Expr::hstack(this->l, Expr::mul(Expr::add(Expr::sub(z3, Expr::mul(3, z2)), Expr::sub(Expr::mul(3, z1), z0)), Qptr)), Domain::inQCone());
		}
		else
		{
			M->constraint(Expr::hstack(Expr::constTerm(nE, 1.), this->l, Expr::mul(Expr::add(Expr::sub(z3, Expr::mul(3, z2)), Expr::sub(Expr::mul(3, z1), z0)), Qptr)), Domain::inRotatedQCone());
		}
	}

	std::shared_ptr<ndarray<int, 1>>  idxsOut = new_array_ptr<int>(this->g->findOutEdges(this->startKey));
	std::shared_ptr<ndarray<int, 1>>  idxtIn = new_array_ptr<int>(this->g->findInEdges(this->targetKey));
	M->constraint(Expr::sum(y->pick(idxsOut)), Domain::equalsTo(1.0));
	M->constraint(Expr::sum(y->pick(idxtIn)), Domain::equalsTo(1.0));

	Variable::t apS, apT, appS, appT, apppS, apppT;
	if ((cpStart->size() > 0) && (this->params.continuityDerivativeOrder >= 1))
	{
		apS = M->variable("apS", 1, Domain::greaterThan(0.0));
	}
	if ((cpTarget->size() > 0) && (this->params.continuityDerivativeOrder >= 1))
	{
		apT = M->variable("apT", 1, Domain::greaterThan(0.0));
	}
	if ((cppStart->size() > 0) && (this->params.continuityDerivativeOrder >= 2))
	{
		appS = M->variable("appS", 1, Domain::greaterThan(0.0));
	}
	if ((cppTarget->size() > 0) && (this->params.continuityDerivativeOrder >= 2))
	{
		appT = M->variable("appT", 1, Domain::greaterThan(0.0));
	}
	if ((cpppStart->size() > 0) && (this->params.continuityDerivativeOrder >= 3))
	{
		apppS = M->variable("apppS", 1, Domain::greaterThan(0.0));
	}
	if ((cpppTarget->size() > 0) && (this->params.continuityDerivativeOrder >= 3))
	{
		apppT = M->variable("apppT", 1, Domain::greaterThan(0.0));
	}

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
			Expression::t exprIn(Expr::zeros(n));
			Expression::t exprInp(Expr::zeros(n));
			Expression::t exprInpp(Expr::zeros(n));
			Expression::t exprInppp(Expr::zeros(n));
			Expression::t exprOut(Expr::zeros(n));
			Expression::t exprOutp(Expr::zeros(n));
			Expression::t exprOutpp(Expr::zeros(n));
			Expression::t exprOutppp(Expr::zeros(n));

			for (int j = 0; j < inEdges.size(); j++)
			{
				//Last control point
				Variable::t z0 = z->slice(new_array_ptr<int, 1>({ (this->N - 1) * n,inEdges[j] }), new_array_ptr<int, 1>({ this->N * n,inEdges[j] + 1 }));
				exprIn = Expr::add(exprIn, z0);
				if (this->params.continuityDerivativeOrder >= 1)
				{
					//First-order derivatives for last control points
					Variable::t z1 = z->slice(new_array_ptr<int, 1>({ (this->N - 2) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 1) * n,inEdges[j] + 1 }));
					exprInp = Expr::add(exprInp, Expr::sub(z0, z1));
				}
				if (this->params.continuityDerivativeOrder >= 2)
				{
					//Second-order derivatives for last control points
					Variable::t z1 = z->slice(new_array_ptr<int, 1>({ (this->N - 2) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 1) * n,inEdges[j] + 1 }));
					Variable::t z2 = z->slice(new_array_ptr<int, 1>({ (this->N - 3) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 2) * n,inEdges[j] + 1 }));
					exprInpp = Expr::add(exprInpp, Expr::add(Expr::sub(z0, Expr::mul(2, z1)), z2));
				}
				if (this->params.continuityDerivativeOrder >= 3)
				{
					//Third-order derivatives for last control points
					Variable::t z1 = z->slice(new_array_ptr<int, 1>({ (this->N - 2) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 1) * n,inEdges[j] + 1 }));
					Variable::t z2 = z->slice(new_array_ptr<int, 1>({ (this->N - 3) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 2) * n,inEdges[j] + 1 }));
					Variable::t z3 = z->slice(new_array_ptr<int, 1>({ (this->N - 4) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 4) * n,inEdges[j] + 1 }));
					exprInppp = Expr::add(exprInppp, Expr::add(Expr::sub(z0, Expr::mul(3, z1)), Expr::sub(Expr::mul(3, z2), z3)));
				}
			}

			for (int j = 0; j < outEdges.size(); j++)
			{
				//First control point
				Variable::t z0 = z->slice(new_array_ptr<int, 1>({ 0,outEdges[j] }), new_array_ptr<int, 1>({ n,outEdges[j] + 1 }));
				exprOut = Expr::add(exprOut, z0);
				if (this->params.continuityDerivativeOrder >= 1)
				{
					//First-order derivatives for first control points
					Variable::t z1 = z->slice(new_array_ptr<int, 1>({ n,outEdges[j] }), new_array_ptr<int, 1>({ 2 * n,outEdges[j] + 1 }));
					exprOutp = Expr::add(exprOutp, Expr::sub(z1, z0));
				}
				if (this->params.continuityDerivativeOrder >= 2)
				{
					//Second-order derivatives for first control points
					Variable::t z1 = z->slice(new_array_ptr<int, 1>({ n,outEdges[j] }), new_array_ptr<int, 1>({ 2 * n,outEdges[j] + 1 }));
					Variable::t z2 = z->slice(new_array_ptr<int, 1>({ 2 * n,outEdges[j] }), new_array_ptr<int, 1>({ 3 * n,outEdges[j] + 1 }));
					exprOutpp = Expr::add(exprOutpp, Expr::add(Expr::sub(z2, Expr::mul(2, z1)), z0));
				}
				if (this->params.continuityDerivativeOrder >= 3)
				{
					//Third-order derivatives for first control points
					Variable::t z1 = z->slice(new_array_ptr<int, 1>({ n,outEdges[j] }), new_array_ptr<int, 1>({ 2 * n,outEdges[j] + 1 }));
					Variable::t z2 = z->slice(new_array_ptr<int, 1>({ 2 * n,outEdges[j] }), new_array_ptr<int, 1>({ 3 * n,outEdges[j] + 1 }));
					Variable::t z3 = z->slice(new_array_ptr<int, 1>({ 3 * n,outEdges[j] }), new_array_ptr<int, 1>({ 4 * n,outEdges[j] + 1 }));
					exprOutppp = Expr::add(exprOutppp, Expr::add(Expr::sub(z3, Expr::mul(3, z2)), Expr::sub(Expr::mul(3, z1), z0)));
				}
			}
			M->constraint(Expr::sub(Expr::sum(y->pick(idxvIn)), Expr::sum(y->pick(idxvOut))), Domain::equalsTo(0.0));
			M->constraint(Expr::sub(exprIn, exprOut), Domain::equalsTo(0.0));
			if (this->params.continuityDerivativeOrder >= 1)
			{
				M->constraint(Expr::sub(exprInp, exprOutp), Domain::equalsTo(0.0));
			}
			if (this->params.continuityDerivativeOrder >= 2)
			{
				M->constraint(Expr::sub(exprInpp, exprOutpp), Domain::equalsTo(0.0));
			}
			if (this->params.continuityDerivativeOrder >= 3)
			{
				M->constraint(Expr::sub(exprInppp, exprOutppp), Domain::equalsTo(0.0));
			}
		}
		else
		{
			if (nodeKeys[i] == startKey)
			{
				std::vector<int> outEdges = this->g->findOutEdges(nodeKeys[i]);
				if (cpStart->size() > 0)
				{
					Expression::t exprSp(Expr::zeros(n));
					if (this->params.continuityDerivativeOrder >= 1)
					{
						//First-order derivatives for first control points
						for (int j = 0; j < outEdges.size(); j++)
						{
							Variable::t z1 = z->slice(new_array_ptr<int, 1>({ n,outEdges[j] }), new_array_ptr<int, 1>({ 2 * n,outEdges[j] + 1 }));
							Variable::t z0 = z->slice(new_array_ptr<int, 1>({ 0,outEdges[j] }), new_array_ptr<int, 1>({ n,outEdges[j] + 1 }));
							exprSp = Expr::add(exprSp, Expr::sub(z1, z0));
						}
						std::cout << Expr::mul(Expr::repeat(apS, 1, 2), cpStart)->getSize() << std::endl;
						//std::cout << Expr::repeat(apS, 1,2)->getSize() << std::endl;
						//std::cout << Expr::mulElm(cpStart,Expr::repeat(apS,2,1))->toString() << std::endl;
						//std::cout << Expr::sub(exprSp, Expr::mul(apS, cpStart))->toString() << std::endl;
						M->constraint(Expr::sub(exprSp, Expr::mul(Expr::repeat(apS, 1, 2), cpStart)), Domain::equalsTo(0.0));
					}
				}
				if (cppStart->size() > 0)
				{
					Expression::t exprSpp(Expr::zeros(n));
					if (this->params.continuityDerivativeOrder >= 2)
					{
						//Second-order derivatives for first control points
						for (int j = 0; j < outEdges.size(); j++)
						{
							Variable::t z2 = z->slice(new_array_ptr<int, 1>({ 2 * n,outEdges[j] }), new_array_ptr<int, 1>({ 3 * n,outEdges[j] + 1 }));
							Variable::t z1 = z->slice(new_array_ptr<int, 1>({ n,outEdges[j] }), new_array_ptr<int, 1>({ 2 * n,outEdges[j] + 1 }));
							Variable::t z0 = z->slice(new_array_ptr<int, 1>({ 0,outEdges[j] }), new_array_ptr<int, 1>({ n,outEdges[j] + 1 }));
							exprSpp = Expr::add(exprSpp, Expr::add(Expr::sub(z2, Expr::mul(2, z1)), z0));
						}
						M->constraint(Expr::sub(exprSpp, Expr::mul(appS, cppStart)), Domain::equalsTo(0.0));
					}
				}
				if (cpppStart->size() > 0)
				{
					Expression::t exprSppp(Expr::zeros(n));
					if (this->params.continuityDerivativeOrder >= 3)
					{
						//Second-order derivatives for first control points
						for (int j = 0; j < outEdges.size(); j++)
						{
							Variable::t z3 = z->slice(new_array_ptr<int, 1>({ 3 * n,outEdges[j] }), new_array_ptr<int, 1>({ 4 * n,outEdges[j] + 1 }));
							Variable::t z2 = z->slice(new_array_ptr<int, 1>({ 2 * n,outEdges[j] }), new_array_ptr<int, 1>({ 3 * n,outEdges[j] + 1 }));
							Variable::t z1 = z->slice(new_array_ptr<int, 1>({ n,outEdges[j] }), new_array_ptr<int, 1>({ 2 * n,outEdges[j] + 1 }));
							Variable::t z0 = z->slice(new_array_ptr<int, 1>({ 0,outEdges[j] }), new_array_ptr<int, 1>({ n,outEdges[j] + 1 }));
							exprSppp = Expr::add(exprSppp, Expr::add(Expr::sub(z3, Expr::mul(3, z2)), Expr::sub(Expr::mul(3, z1), z0)));
						}
						M->constraint(Expr::sub(exprSppp, Expr::mul(apppS, cpppStart)), Domain::equalsTo(0.0));
					}
				}
			}
			if (nodeKeys[i] == targetKey)
			{
				std::vector<int> inEdges = this->g->findInEdges(nodeKeys[i]);
				if (cpTarget->size() > 0)
				{
					Expression::t exprTp(Expr::zeros(n));
					if (this->params.continuityDerivativeOrder >= 1)
					{
						//First-order derivatives for last control points
						for (int j = 0; j < inEdges.size(); j++)
						{
							Variable::t z1 = z->slice(new_array_ptr<int, 1>({ (this->N - 1) * n,inEdges[j] }), new_array_ptr<int, 1>({ this->N * n,inEdges[j] + 1 }));
							Variable::t z0 = z->slice(new_array_ptr<int, 1>({ (this->N - 2) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 1) * n,inEdges[j] + 1 }));
							exprTp = Expr::add(exprTp, Expr::sub(z1, z0));
						}
						M->constraint(Expr::sub(exprTp, Expr::mul(apT, cpTarget)), Domain::equalsTo(0.0));
					}
				}
				if (cppTarget->size() > 0)
				{
					Expression::t exprTpp(Expr::zeros(n));
					if (this->params.continuityDerivativeOrder >= 1)
					{
						//Second-order derivatives for last control points
						for (int j = 0; j < inEdges.size(); j++)
						{
							Variable::t z2 = z->slice(new_array_ptr<int, 1>({ (this->N - 1) * n,inEdges[j] }), new_array_ptr<int, 1>({ this->N * n,inEdges[j] + 1 }));
							Variable::t z1 = z->slice(new_array_ptr<int, 1>({ (this->N - 2) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 1) * n,inEdges[j] + 1 }));
							Variable::t z0 = z->slice(new_array_ptr<int, 1>({ (this->N - 3) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 2) * n,inEdges[j] + 1 }));
							exprTpp = Expr::add(exprTpp, Expr::add(Expr::sub(z2, Expr::mul(2, z1)), z0));
						}
						M->constraint(Expr::sub(exprTpp, Expr::mul(appT, cppTarget)), Domain::equalsTo(0.0));
					}
				}
				if (cpppTarget->size() > 0)
				{
					Expression::t exprTppp(Expr::zeros(n));
					if (this->params.continuityDerivativeOrder >= 1)
					{
						//Third-order derivatives for last control points
						for (int j = 0; j < inEdges.size(); j++)
						{
							Variable::t z3 = z->slice(new_array_ptr<int, 1>({ (this->N - 1) * n,inEdges[j] }), new_array_ptr<int, 1>({ this->N * n,inEdges[j] + 1 }));
							Variable::t z2 = z->slice(new_array_ptr<int, 1>({ (this->N - 2) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 1) * n,inEdges[j] + 1 }));
							Variable::t z1 = z->slice(new_array_ptr<int, 1>({ (this->N - 3) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 2) * n,inEdges[j] + 1 }));
							Variable::t z0 = z->slice(new_array_ptr<int, 1>({ (this->N - 4) * n,inEdges[j] }), new_array_ptr<int, 1>({ (this->N - 3) * n,inEdges[j] + 1 }));
							exprTppp = Expr::add(exprTppp, Expr::add(Expr::sub(z3, Expr::mul(3, z2)), Expr::sub(Expr::mul(3, z1), z0)));
						}
						M->constraint(Expr::sub(exprTppp, Expr::mul(apppT, cpppTarget)), Domain::equalsTo(0.0));
					}
				}
			}
		}
	}
	std::vector<Edge> edges(this->g->getEdges());
	int i = 0;
	for (std::vector<Edge>::iterator edgeIt = edges.begin(); edgeIt != edges.end(); edgeIt++)
	{
		int u = edgeIt->first;
		int v = edgeIt->second;

		//First control point
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

		//Control points in the middle
		std::map<std::pair<int, int>, int>::iterator itNavGraph2GCS = this->navGraph2gcs->find(*edgeIt);
		if (itNavGraph2GCS != this->navGraph2gcs->end())
		{
			int gcsKey = itNavGraph2GCS->second;
			PolyhedronNode* node = (PolyhedronNode*)this->gcs->getNode(gcsKey);
			std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(node->polyhedron.A);
			std::shared_ptr<ndarray<double, 2>> b_ptr = Eigen2NdArray(Eigen::MatrixXd(node->polyhedron.b));
			for (int j = 1; j < (N - 1); j++)
			{
				Expression::t z_expr = z->slice(new_array_ptr<int, 1>({ j * n,i }), new_array_ptr<int, 1>({ (j + 1) * n,i + 1 }));
				M->constraint(Expr::sub(Expr::mul(A_ptr, z_expr), Expr::mul(b_ptr, y->pick(new_array_ptr<int, 1>({ i })))), Domain::lessThan(0.0));
			}
		}
		else
		{
			Range* range = Range::getInstance();
			Eigen::MatrixXd ARange;
			Eigen::VectorXd bRange;
			range->getConstraints(ARange, bRange);
			std::shared_ptr<ndarray<double, 2>> A_ptr = Eigen2NdArray(ARange);
			std::shared_ptr<ndarray<double, 2>> b_ptr = Eigen2NdArray(Eigen::MatrixXd(bRange));
			for (int j = 1; j < (N - 1); j++)
			{
				Expression::t z_expr = z->slice(new_array_ptr<int, 1>({ j * n,i }), new_array_ptr<int, 1>({ (j + 1) * n,i + 1 }));
				M->constraint(Expr::sub(Expr::mul(A_ptr, z_expr), Expr::mul(b_ptr, y->pick(new_array_ptr<int, 1>({ i })))), Domain::lessThan(0.0));
			}
		}
		//Last control point
		if (v == targetKey)
		{
			Expression::t z_expr = z->slice(new_array_ptr<int, 1>({ (this->N - 1) * n,i }), new_array_ptr<int, 1>({ this->N * n,i + 1 }));
			Expression::t qy = Expr::mul(this->qtarget, y->pick(new_array_ptr<int, 1>({ i })));
			M->constraint(Expr::sub(z_expr, qy), Domain::equalsTo(0.0));
		}
		else
		{
			PolyhedronNode* vNode = (PolyhedronNode*)this->g->getNode(v);
			Expression::t z_expr = z->slice(new_array_ptr<int, 1>({ (this->N - 1) * n,i }), new_array_ptr<int, 1>({ this->N * n,i + 1 }));
			std::shared_ptr<ndarray<double, 2>> Ap_ptr = Eigen2NdArray(vNode->polyhedron.A);
			std::shared_ptr<ndarray<double, 2>> bp_ptr = Eigen2NdArray(Eigen::MatrixXd(vNode->polyhedron.b));
			M->constraint(Expr::sub(Expr::mul(Ap_ptr, z_expr), Expr::mul(bp_ptr, y->pick(new_array_ptr<int, 1>({ i })))), Domain::lessThan(0.0));
		}
		i++;
	}
}

void ConvexRelaxationBezierCurveSPP_GCS::setObjective(Model::t& M)
{
	//M->objective(ObjectiveSense::Minimize, Expr::sum(this->l));
	M->objective(ObjectiveSense::Minimize, Expr::add(Expr::sum(this->l), Expr::mul(Expr::sum(this->y), 1e-4)));
}

void ConvexRelaxationBezierCurveSPP_GCS::optimize()
{
	Model::t M = new Model("ConvexRelaxationBezierCurve");
	auto _M = finally([&]() { M->dispose(); });
	int nE = this->g->numEdges;
	int nN = this->g->numNodes;
	this->allocateVariables(M);
	this->setConstraints(M);
	this->setObjective(M);
	//M->writeTask("dump.ptf");

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
}

Eigen::MatrixXd ConvexRelaxationBezierCurveSPP_GCS::computeBezierL2NormMatrix()
{
	int NN = this->N - this->params.costDerivativeOrder;
	assert(NN >= 2); //It does not make sense to call this function with such few control points
	Eigen::MatrixXd Q(NN, NN);
	Eigen::VectorXd C1(NN);
	Eigen::VectorXd C2(2 * NN - 1);
	for (int n = 0; n <= (NN - 1); n++)
		C1(n) = nchoosek(NN - 1, n);
	for (int n = 0; n <= 2 * (NN - 1); n++)
		C2(n) = nchoosek(2 * (NN - 1), n);
	for (int n = 0; n <= (NN - 1); n++)
	{
		for (int m = 0; m <= (NN - 1); m++)
		{
			Q(n, m) = C1(n) * C1(m) / (C2(n + m) * (2 * N - 1));
		}
	}
	Eigen::MatrixXd Qbig(NN * this->n, NN * this->n);
	Eigen::MatrixXd I = Eigen::MatrixXd::Identity(this->n, this->n);
	for (int i = 0; i < NN; i++)
	{
		for (int j = 0; j < NN; j++)
		{
			Qbig.block(i * this->n, j * this->n, this->n, this->n) = I * Q(i, j);
		}
	}
	if (this->params.squaredRootL2NormCost)
	{
		Qbig.llt();
		//Now update the pointer
		this->Qptr = Eigen2NdArray(Qbig);
	}
	else
	{
		this->Qptr = Eigen2NdArray(Qbig);
	}
	return Qbig;
}

ConvexRelaxationBezierCurveParams_t ConvexRelaxationBezierCurveSPP_GCS::getDefaultBezierCurveSolverParams()
{
	ConvexRelaxationBezierCurveParams_t params = {5,1,1,false,false};
	return params;
}