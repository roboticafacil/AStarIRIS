#include "CObsConic.h"

CObsConic::CObsConic()
{
}

CObsConic::CObsConic(CObsConic* cObs)
{
	this->conicObjects = cObs->conicObjects;
	this->numObjects = cObs->numObjects;
}

void CObsConic::addObject(ConicSet* object)
{
	conicObjects.push_back(object);
	numObjects++;
}

bool CObsConic::isFree(const Eigen::VectorXd& q, const double &tol)
{
	for (std::vector<ConicSet*>::iterator it = conicObjects.begin(); it != conicObjects.end(); it++)
	{
		ConicSet* conicSet = *it;
		if (conicSet->isInside(q,tol))
			return false;
	}
	return true;
}

void CObsConic::print(std::ostream& out)
{
	for (int i = 0; i < this->conicObjects.size(); i++)
	{
		ConicSet* conicSet = this->conicObjects[i];
		if (dynamic_cast<Polyhedron*>(conicSet) != NULL)
		{
			Polyhedron* poly = (Polyhedron*)conicSet;
			poly->print(out, "A", "b",true);
			out << "AObs{" << i + 1 << "}=-A;" << std::endl;
			out << "bObs{" << i + 1 << "}=-b;" << std::endl;
			out << "colObs{" << i + 1 << "}='r';" << std::endl;
		}
	}
}
