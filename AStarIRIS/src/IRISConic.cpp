#include "IRISConic.h"
#include "Ellipsoid.h"
#include "PolyhedronNode.h"
#include "IRISConic.h"


IRISConic::IRISConic(CObsConic& cObs, const IRISParams_t& params): _cObs(&cObs), params(params)
{	
}

void IRISConic::generateGCS()
{
	int trials = 0;
	while (trials < params.maxTrials)
	{
		//if (gcs->numNodes >= 26)
		//	std::cout << "Start debugging" << std::endl;
		Eigen::VectorXd seed = this->generateRandomSeed(trials);
		//std::cout << "Seed" << std::endl;
		//std::cout << seed << std::endl;
		if (trials >= params.maxTrials)
			break;
		
		addConvexSets(seed);  //When we call this function, we allocate memory of the generated convexSets. Where we should delete them?
	}
}

Eigen::VectorXd IRISConic::generateRandomSeed(int& trials)
{
	Range* range = Range::getInstance();
	std::pair<double, double> bounds = range->getBounds();
	Eigen::VectorXd seed;
	while (true)
	{
		seed = Eigen::VectorXd::Random(params.n, 1);
		seed = (seed + Eigen::VectorXd::Constant(params.n, 1, 1.)) * (bounds.second - bounds.first) / 2.;
		seed = (seed + Eigen::VectorXd::Constant(params.n, 1, bounds.first));
		if (!this->gcs.contains(seed,params.tol))
		{
			if (this->_cObs->isFree(seed,params.tol))
			{
				trials = 0;
				return seed;
			}
		}
		trials++;
		if (trials >= params.maxTrials)
			break;
	}
}

void IRISConic::addConvexSets(const Eigen::VectorXd& q)
{
	//Warning, convexSet is not removed from memory!!
	Ellipsoid ellipsoid(Eigen::MatrixXd::Identity(q.rows(), q.rows()), q);
	Polyhedron* convexSet = new Polyhedron(q.rows());
	std::vector<IRISNeighbour_t> neighbours;
	while(true)
	{
		this->computeConvexSet(ellipsoid, *convexSet, neighbours);
		//Add the convexSet to GCS
		Eigen::VectorXd bShrinked = (1.-params.shrinkFactor)*convexSet->b + params.shrinkFactor * (convexSet->A*ellipsoid.getCentroid());
		int nodeKey=this->gcs.addNode(new PolyhedronNode(convexSet->A, convexSet->b,bShrinked));
		//Check if the original seed is inside the generated convex set. If not, repeat
		if (convexSet->isInside(q))
			break;
		ellipsoid = Ellipsoid(Eigen::MatrixXd::Identity(q.rows(), q.rows()), q);
	}
	//convexSet->print();
}

void IRISConic::computeConvexSet(Ellipsoid& ellipsoid, Polyhedron& convexSet, std::vector<IRISNeighbour_t>& neighbours)
{
	while (true)
	{
		separatingHyperplanes(ellipsoid, convexSet, neighbours);
		//convexSet->print();
		Ellipsoid newEllipsoid = convexSet.inscribedEllipsoid();
		//newEllipsoid.print();
		double detC = ellipsoid.C.determinant();
		double detNewC = newEllipsoid.C.determinant();
		if (((detNewC - detC) / detC) < params.tolConvexSetConvergence)
		{
			ellipsoid = newEllipsoid;
			break;
		}
		ellipsoid = newEllipsoid;
	}
	std::vector<int> removedConstraints=convexSet.removeRepeatedConstraints();
	for (std::vector<int>::iterator it1 = removedConstraints.begin(); it1 != removedConstraints.end(); it1++)
	{
		for (std::vector<IRISNeighbour_t>::iterator it = neighbours.begin(); it != neighbours.end();)
		{
			if (it->idx == *it1)
				neighbours.erase(it);
			else
				it++;
		}
	}
}

void IRISConic::separatingHyperplanes(Ellipsoid& ellipsoid, Polyhedron& polyhedron, std::vector<IRISNeighbour_t>& neighbours)
{
	neighbours.clear();
	Eigen::VectorXd centroid = ellipsoid.getCentroid();
	int n = centroid.rows();
	std::vector<Node*> nodes = gcs.getNodes();
	std::vector<int> nodeKeys = gcs.getNodeKeys();
	std::vector<ConicSet*> objects;
	std::vector<int> objectsIdx;
	std::vector<double> distances;
	std::vector<Eigen::VectorXd> closestPoints;
	std::vector<int> graphNodeKeys;
	int idx = 0;
	for (std::vector<ConicSet*>::iterator it = _cObs->conicObjects.begin(); it != _cObs->conicObjects.end(); it++)
	{
		ConicSet* p = *it;
		Eigen::VectorXd p_out(n);
		double d = p->closestPointExpandingEllipsoid(ellipsoid,p_out);
		objects.push_back(p);
		objectsIdx.push_back(idx++);
		distances.push_back(d);
		closestPoints.push_back(p_out);
		graphNodeKeys.push_back(-1);
	}
	int i = 0;
	for (std::vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		Node* node = *it;
		PolyhedronNode* poly = (PolyhedronNode*)node->getNodeData();
		Eigen::VectorXd p_out(n);
		double d;
		if (poly->useShrinked)
		{
			d = poly->shrinkedPolyhedron.closestPointExpandingEllipsoid(ellipsoid, p_out);
			objects.push_back(&poly->shrinkedPolyhedron);
		}
		else
		{
			d = poly->polyhedron.closestPointExpandingEllipsoid(ellipsoid, p_out);
			objects.push_back(&poly->polyhedron);
		}
		objectsIdx.push_back(idx++);
		distances.push_back(d);
		closestPoints.push_back(p_out);
		graphNodeKeys.push_back(nodeKeys[i]);
		i++;
	}
	Eigen::MatrixXd A(0,n);
	Eigen::VectorXd b(0);
	while (!objects.empty())
	{
		std::vector<double>::iterator min_dist_it = std::min_element(distances.begin(), distances.end());
		double min_dist=*min_dist_it; //Min distance
		int ii = std::distance(distances.begin(),min_dist_it); //Index of the element with the minimum distance
		//std::cout << "Closest object " << objectsIdx[ii] << std::endl;
		//std::cout << "Closest object distance " << min_dist << std::endl;
		Eigen::VectorXd x = closestPoints[ii];
		//std::cout << "Closest object point" << std::endl;
		//std::cout << x << std::endl;
		Eigen::VectorXd ai(n);
		double bi;
		ellipsoid.tangentHyperplane(x, ai, bi);
		//std::cout << "Closest object separting hyperplane " << std::endl;
		//std::cout << "ai=";
		//std::cout << ai << std::endl;
		//std::cout << "bi=";
		//std::cout << bi << std::endl;
		//Append the hyperplane to the inequalities constraints
		A.conservativeResize(A.rows() + 1, A.cols());
		A.row(A.rows()-1) = ai.transpose();
		b.conservativeResize(b.rows() + 1);
		b(b.rows() - 1) = bi;
		if (graphNodeKeys[ii] >= 0)
		{
			//std::cout << "Object added to the neighbours list with key" << graphNodeKeys[ii] << std::endl;
			IRISNeighbour_t neighbour;
			neighbour.nodeKey = graphNodeKeys[ii];
			neighbour.ai = ai;
			neighbour.bi = bi;
			neighbour.idx = A.rows() - 1;
			neighbours.push_back(neighbour);
		}

		//Delete the object from the lists
		objects.erase(objects.begin() + ii);
		objectsIdx.erase(objectsIdx.begin() + ii);
		distances.erase(distances.begin() + ii);
		closestPoints.erase(closestPoints.begin() + ii);
		graphNodeKeys.erase(graphNodeKeys.begin() + ii);
		//Remove objects that are outside the separating hyperplane
		std::vector<int>::iterator itIdx = objectsIdx.begin();
		for (std::vector<ConicSet*>::iterator it=objects.begin();it!=objects.end();)
		{
			ConicSet* c = *it;
			bool isInside = c->isInsideSeparatingHyperplane(ai, bi,params.tol);
			//std::cout << "Object " << *itIdx << "is inside?" << isInside << std::endl;
			int j = std::distance(objects.begin(), it); 
			if (!isInside)
			{
				//std::cout << "Original object removed by separating hyperplane " << *itIdx << std::endl;
				objects.erase(objects.begin() + j);
				objectsIdx.erase(objectsIdx.begin() + j);
				distances.erase(distances.begin() + j);
				closestPoints.erase(closestPoints.begin() + j);
				graphNodeKeys.erase(graphNodeKeys.begin() + j);
			}
			else
			{
				it++; //Increase the iterator if not removed
				itIdx++;
			}
			//std::cout << it._Ptr << std::endl;
			//std::cout << objects.end()._Ptr << std::endl;
			if (it._Ptr == objects.end()._Ptr)
				break;
			
		}
	}
	polyhedron.update(A,b);
}

IRISParams_t IRISConic::getDefaultIRISParams()
{
	IRISParams_t params = {2,1e-2,0.25,999,1e-3};
	return params;
}

/*void IRISConic::addStartNode(const Eigen::VectorXd& q)
{
	this->addConvexSets(q);
	int startKeyGCS = this->gcs.findConvexSet(q);
	PolyhedronNode* startNodeGCS = (PolyhedronNode*)this->gcs.getNode(startKeyGCS)->getNodeData();
	for (int i = 0; i < startNodeGCS->polyhedron.A.rows(); i++)
	{
		Eigen::VectorXd ai = Eigen::VectorXd(startNodeGCS->polyhedron.A.row(i));
		double bi = startNodeGCS->polyhedron.b(i);
	}
}*/

