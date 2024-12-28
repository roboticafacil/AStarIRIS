#include "IRISConic.h"
#include "Ellipsoid.h"
#include "PolyhedronNode.h"
#include "IRISConic.h"


IRISConic::IRISConic(CObsConic* cObs, GCS* gcs, const IRISParams_t& params): gcs(gcs), cObs(cObs), params(params)
{	
}

void IRISConic::generateGCS()
{
	int count = 0;
	while (count < params.maxCount)
	{
		Eigen::VectorXd seed = this->generateRandomSeed(count);
		//std::cout << "Seed" << std::endl;
		//std::cout << seed << std::endl;
		if (count >= params.maxCount)
			break;
		addConvexSets(seed);
	}
}

Eigen::VectorXd IRISConic::generateRandomSeed(int& count)
{
	Range* range = Range::getInstance();
	std::pair<double, double> bounds = range->getBounds();
	Eigen::VectorXd seed;
	while (true)
	{
		seed = Eigen::VectorXd::Random(params.n, 1);
		seed = (seed + Eigen::VectorXd::Constant(params.n, 1, 1.)) * (bounds.second - bounds.first) / 2.;
		seed = (seed + Eigen::VectorXd::Constant(params.n, 1, bounds.first));
		if (!gcs->contains(seed))
		{
			if (cObs->isFree(seed))
			{
				count = 0;
				return seed;
			}
		}
		count++;
		if (count >= params.maxCount)
			break;
	}
}

void IRISConic::addConvexSets(const Eigen::VectorXd& q)
{
	Ellipsoid ellipsoid(Eigen::MatrixXd::Identity(q.rows(), q.rows()), q);
	Polyhedron* convexSet = new Polyhedron(q.rows());
	std::vector<IRISNeighbour_t> neighbours;
	std::vector<int> newNodeKeys;
	//convexSet.allocateClosestPointEllipsoidSolver();
	//convexSet.allocateInscribedEllipsoidSolver();
	while(true)
	{
		while (true)
		{
			separatingHyperplanes(ellipsoid,convexSet, neighbours);
			//convexSet->print();
			Ellipsoid newEllipsoid = convexSet->inscribedEllipsoid();
			//newEllipsoid.print();
			double detC = ellipsoid.C.determinant();
			double detNewC = newEllipsoid.C.determinant();
			if (((detNewC - detC) / detC) < params.tolConvexSetConvergence)
				break;
			ellipsoid = newEllipsoid;
		}
		Eigen::VectorXd bShrinked = (1.-params.shrinkFactor)*convexSet->b + params.shrinkFactor * (convexSet->A*ellipsoid.getCentroid());
		int nodeKey=gcs->addNode(new PolyhedronNode(convexSet->A, convexSet->b,bShrinked));
		newNodeKeys.push_back(nodeKey);
		for (std::vector<IRISNeighbour_t>::iterator it = neighbours.begin(); it != neighbours.end(); it++)
		{
			gcs->addEdge(nodeKey,it->nodeKey);
			gcs->addEdge(it->nodeKey,nodeKey);
		}
		if (convexSet->isInside(q))
			break;
		ellipsoid = Ellipsoid(Eigen::MatrixXd::Identity(q.rows(), q.rows()), q);
	}
	//convexSet->print();
	//Now they have generated, we can treat them as obstacles for nexts calls
	/*for (std::vector<int>::iterator it = newNodeKeys.begin(); it != newNodeKeys.end(); it++)
	{
		Node* node = gcs->getNode(*it);
		PolyhedronNode* polyNode=(PolyhedronNode*)node->getNodeData();
		polyNode->useShrinked = false;
	}*/
}

void IRISConic::separatingHyperplanes(Ellipsoid& ellipsoid, Polyhedron* polyhedron, std::vector<IRISNeighbour_t>& neighbours)
{
	neighbours.clear();
	Eigen::VectorXd centroid = ellipsoid.getCentroid();
	int n = centroid.rows();
	std::vector<Node*> nodes = gcs->getNodes();
	std::vector<int> nodeKeys = gcs->getNodeKeys();
	std::vector<ConicSet*> objects;
	std::vector<int> objectsIdx;
	std::vector<double> distances;
	std::vector<Eigen::VectorXd> closestPoints;
	std::vector<int> graphNodeKeys;
	int idx = 0;
	for (std::vector<ConicSet*>::iterator it = cObs->conicObjects.begin(); it != cObs->conicObjects.end(); it++)
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
			if (!c->isInsideSeparatingHyperplane(ai, bi))
			{
				//std::cout << "Original object removed by separating hyperplane " << *itIdx << std::endl;
				int j = std::distance(objects.begin(), it);
				objects.erase(objects.begin() + j);
				objectsIdx.erase(objectsIdx.begin() + j);
				distances.erase(distances.begin() + j);
				closestPoints.erase(closestPoints.begin() + j);
				graphNodeKeys.erase(graphNodeKeys.begin() + j);
			}
			else
				it++; //Increase the iterator if not removed
			//std::cout << it._Ptr << std::endl;
			//std::cout << objects.end()._Ptr << std::endl;
			if (it._Ptr == objects.end()._Ptr)
				break;
			itIdx++;
		}
	}
	polyhedron->update(A,b);
}

IRISParams_t IRISConic::getDefaultParams()
{
	IRISParams_t params = {2,1e-2,0.25,99};
	return params;
}

