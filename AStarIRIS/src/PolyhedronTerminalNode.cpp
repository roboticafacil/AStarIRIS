#include<Eigen/Dense>
#include "Polyhedron.h"
#include "Node.h"
#include "Range.h"
#include "PolyhedronTerminalNode.h"

PolyhedronTerminalNode::PolyhedronTerminalNode(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const int& terminalConstrainIdx) : PolyhedronNode(A,b), terminalConstrainIdx(terminalConstrainIdx), expanded(false)
{
	this->n = this->polyhedron.A.cols();
	this->m = this->polyhedron.A.rows();
	Eigen::VectorXd ai = this->polyhedron.A.row(this->terminalConstrainIdx);
	this->bi = this->polyhedron.b(this->terminalConstrainIdx);
	this->aipivot=ai.cwiseAbs().maxCoeff(&this->col);
	for (int i = 0; i < this->m; i++)
	{
		if (i!=terminalConstrainIdx)
			this->rows.push_back(i);
	}
	for (int i = 0; i < this->n; i++)
	{
		if (i != col)
			this->cols.push_back(i);
	}
	this->Aj = this->polyhedron.A(rows,Eigen::placeholders::all);
	this->bj = this->polyhedron.b(rows);
	this->ai_col = ai(cols);
}

void* PolyhedronTerminalNode::getNodeData() {
	return this;
}

bool PolyhedronTerminalNode::generateRandomSeed(Eigen::VectorXd& seed, const double &tol, const int& maxTrials)
{
	bool contained_Ajbj = false;
	Range* range = Range::getInstance();
	std::pair<double, double> bounds = range->getBounds();
	Eigen::MatrixXd Arange;
	Eigen::VectorXd brange;
	range->getConstraints(Arange, brange);
	int trials = 0;
	Ellipsoid ellipsoid = this->polyhedron.inscribedEllipsoid();
	while (!contained_Ajbj)
	{
		seed = Eigen::VectorXd::Random(n, 1);
		//Generate a random seend inside the inscribed ellipsoid
		seed.normalize();
		seed=ellipsoid.C*seed+ellipsoid.getCentroid();
		//seed = (seed + Eigen::VectorXd::Constant(n, 1, 1.)) * (bounds.second - bounds.first) / 2.;
		//seed = (seed + Eigen::VectorXd::Constant(n, 1, bounds.first));
		//Project the random seed onto the constraint
		//std::cout << "A=" << this->polyhedron.A << std::endl;
		//std::cout << "b=" << this->polyhedron.b << std::endl;
		//std::cout << "terminalConstraintIdx=" << this->terminalConstrainIdx << std::endl;
		//std::cout << "aipivot=" << this->aipivot << std::endl;
		//std::cout << "ai_col=" << this->ai_col << std::endl;
		//std::cout << "Aj=" << this->Aj << std::endl;
		//std::cout << "bj=" << this->bj << std::endl;
		//std::cout << "seed= " << seed << std::endl;
		//seed(col) = (this->bi+(2.*tol)-this->ai_col.dot(seed(cols))) / this->aipivot;
		//seed(col) = (this->bi -(tol/0.5) - this->ai_col.dot(seed(cols))) / this->aipivot;
		bool cond1 = ((Aj * seed - bj).array() <= 0.).all();
		bool cond2 = ((Arange * seed - brange).array() <= 0.).all();
		//std::cout << "projected seed= " << seed << std::endl;
		//std::cout << "cond1= " << cond1 << std::endl;
		//std::cout << "cond2= " << cond2 << std::endl;
		if ((cond1) && (cond2))
		{
			contained_Ajbj = true;
		}
		trials++;
		if (trials >= maxTrials)
		{
			//std::cout << "Warning: maximum number of trials reached in PolyhedronTerminalNode::generateRandomSeed() " << std::endl;
			break;
		}
	}
	return contained_Ajbj;
}