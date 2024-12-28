#include "Range.h"
#include "fusion.h"
#include "EigenNdArray.h"

using namespace monty;

Range* Range::_range = nullptr;

void Range::setRange(const int &n, const double& lb, const double& ub)
{
	this->_lb = lb;
	this->_ub = ub;
	A.resize(2*n,n);
	b.resize(2*n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
				A(i,j) = 1.;
			else
				A(i,j) = 0.;
		}
		b(i) = ub;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
				A(i+n,j) = -1.;
			else
				A(i+n,j) = 0.;
		}
		b(i+n) = -lb;
	}
	A_ptr = Eigen2NdArray(this->A);
	b_ptr = Eigen2NdArray(this->b);
}

std::pair<double, double> Range::getBounds()
{
	return std::make_pair(this->_lb, this->_ub);
}

std::pair< std::shared_ptr<ndarray<double, 2>>, std::shared_ptr<ndarray<double, 1>> > Range::getConstraints()
{
	return std::make_pair(this->A_ptr, this->b_ptr);
}

Range::Range(): _lb(-1.0), _ub(1.0) //, Polyhedron(0)
{
}

Range* Range::getInstance()
{
	if (_range == nullptr) {
		_range = new Range();
	}
	return _range;
}