#include "Polyhedron.h"
#include "Range.h"
#include "GCS.h"
#include "NavGraph.h"
#include "PolyhedronNode.h"
#include "PointNode.h"
#include "ConvexRelaxationMinDistanceSPP_GCS.h"
#include "MIPMinDistanceSPP_GCS.h"
#include <fstream>

int analyse_convex_relaxation_ambiguities()
{
	std::cout << "Analyse convex relaxation ambiguities" << std::endl;
	Range* range = Range::getInstance();
	double lb = -2.;
	double ub = 2.;
	range->setRange(2, lb, ub);

	std::ofstream fout("analyse_convex_relaxation_ambiguities.m");
	fout << "close all;" << std::endl;
	fout << "clear;" << std::endl;
	fout << "xlim = [" << lb << " " << ub << "];" << std::endl;
	fout << "ylim = [" << lb << " " << ub << "];" << std::endl;
	fout << "Range = PolyHedronRange([xlim; ylim]');" << std::endl;

	GCS gcs = GCS();
	NavGraph navGraph = NavGraph();
	//std::map<std::pair<int, int>, int> navGraph2gcs;
	//std::map<std::pair<int, int>, int> gcs2navGraph;

	PointNode startNode = PointNode(Eigen::Vector<double,2>({-0.5,0.5}));
	int startKey = navGraph.addNode(&startNode);

	PointNode targetNode = PointNode(Eigen::Vector<double, 2>({0.5,-0.5}));
	int targetKey = navGraph.addNode(&targetNode);

	//Let's add the first convex set
	Eigen::Matrix<double,4,2> A1({ {1.,0.},{0.,1.},{-1.,0.},{0.,-1.} });
	Eigen::Vector<double,4> b1({ 0.1,1.1,1.1,0.1 });
	PolyhedronNode c1(A1, b1);
	int c1Key = gcs.addNode(&c1);

	//Now, let's add the second convex set
	Eigen::Matrix<double, 4, 2> A2({ {1.,0.},{0.,1.},{-1.,0.},{0.,-1.} });
	Eigen::Vector<double, 4> b2({1.1,1.1,0.1,0.1 });
	PolyhedronNode c2(A2, b2);
	int c2Key = gcs.addNode(&c2);
	//Add the edge to GCS graph
	gcs.addEdge(c2Key, c1Key);
	//Compute the intersection between c2 and c1
	Polyhedron pI1 = Polyhedron::intersection(Polyhedron(A2, b2),Polyhedron(A1, b1));
	PolyhedronNode I1 = PolyhedronNode(pI1.A, pI1.b);
	//Add the intersection polyhedron to the navGraph
	int I1Key = navGraph.addNode(&I1);
	navGraph.addEdge(startKey, I1Key);
	//navGraph2gcs.emplace(std::make_pair(startKey,I1Key),c1Key);

	//Register the relation betwen the c1-c2 and the intersection in this list
	//gcs2navGraph.emplace(std::make_pair(c2Key,c1Key),I1Key);

	//Now let's add c3
	Eigen::Matrix<double, 4, 2> A3({ {1.,0.},{0.,1.},{-1.,0.},{0.,-1.} });
	Eigen::Vector<double, 4> b3({ 0.1,0.1,1.1,1.1 });
	PolyhedronNode c3(A3, b3);
	int c3Key = gcs.addNode(&c3);
	//Add the edges to GCS graph
	gcs.addEdge(c3Key, c1Key);
	gcs.addEdge(c3Key, c2Key);
	//Compute the intersection between c3 and c1
	Polyhedron pI2 = Polyhedron::intersection(Polyhedron(A3, b3), Polyhedron(A1, b1));
	PolyhedronNode I2 = PolyhedronNode(pI2.A, pI2.b);
	int I2Key = navGraph.addNode(&I2);
	//gcs2navGraph.emplace(std::make_pair(c3Key, c1Key), I2Key);
	navGraph.addEdge(I2Key, I1Key);
	//navGraph2gcs.emplace(std::make_pair(I2Key,I1Key), c1Key);
	navGraph.addEdge(I1Key, I2Key);
	//navGraph2gcs.emplace(std::make_pair(I1Key, I2Key), c1Key);
	navGraph.addEdge(startKey, I2Key);
	//navGraph2gcs.emplace(std::make_pair(startKey, I2Key), c1Key);
	//Compute the intersection between c3 and c2
	Polyhedron pI3 = Polyhedron::intersection(Polyhedron(A3, b3),Polyhedron(A2, b2));
	PolyhedronNode I3 = PolyhedronNode(pI3.A, pI3.b);
	int I3Key = navGraph.addNode(&I3);
	//gcs2navGraph.emplace(std::make_pair(c3Key, c2Key), I3Key);
	navGraph.addEdge(I3Key, I1Key);
	//navGraph2gcs.emplace(std::make_pair(I3Key, I1Key), c1Key);
	navGraph.addEdge(I1Key, I3Key);
	//navGraph2gcs.emplace(std::make_pair(I1Key, I3Key), c1Key);
	navGraph.addEdge(I3Key, I2Key);
	//navGraph2gcs.emplace(std::make_pair(I3Key, I2Key), c2Key);
	navGraph.addEdge(I2Key, I3Key);
	//navGraph2gcs.emplace(std::make_pair(I2Key, I3Key), c2Key);
	navGraph.addEdge(startKey, I3Key);

	//Now let's add c4
	Eigen::Matrix<double, 4, 2> A4({ {1.,0.},{0.,1.},{-1.,0.},{0.,-1.} });
	Eigen::Vector<double, 4> b4({ 1.1,0.1,0.1,1.1 });
	PolyhedronNode c4(A4, b4);
	int c4Key = gcs.addNode(&c4);
	//Add the edges to GCS graph
	gcs.addEdge(c1Key, c4Key);
	gcs.addEdge(c2Key, c4Key);
	gcs.addEdge(c3Key, c4Key);
	navGraph.addEdge(I3Key, targetKey);
	//navGraph2gcs.emplace(std::make_pair(I4Key, targetKey),c4Key);
	//Compute the intersection between c4 and c2
	Polyhedron pI4 = Polyhedron::intersection(Polyhedron(A4, b4), Polyhedron(A2, b2));
	PolyhedronNode I4 = PolyhedronNode(pI4.A, pI4.b);
	int I4Key = navGraph.addNode(&I4);
	//gcs2navGraph.emplace(std::make_pair(c4Key, c2Key), I5Key); //Revisar!!
	navGraph.addEdge(I4Key, I1Key);
	//navGraph2gcs.emplace(std::make_pair(I5Key, I1Key), c2Key); //Revisar!!
	navGraph.addEdge(I1Key, I4Key);
	//gcs2navGraph.emplace(std::make_pair(c4Key, c2Key), I5Key);
	navGraph.addEdge(I4Key, I3Key); 
	//navGraph2gcs.emplace(std::make_pair(I5Key, I3Key), c2Key); //Revisar!!
	navGraph.addEdge(I3Key, I4Key);
	//navGraph2gcs.emplace(std::make_pair(I3Key, I5Key), c2Key); //Revisar!!
	//navGraph2gcs.emplace(std::make_pair(I4Key, I5Key), c2Key); //Revisar!!
	navGraph.addEdge(I4Key, targetKey);
	//navGraph2gcs.emplace(std::make_pair(I5Key, targetKey), c4Key);
	//Compute the intersection between c4 and c3
	Polyhedron pI6 = Polyhedron::intersection(Polyhedron(A4, b4), Polyhedron(A3, b3));
	PolyhedronNode I6 = PolyhedronNode(pI6.A, pI6.b);
	int I6Key = navGraph.addNode(&I6);
	//gcs2navGraph.emplace(std::make_pair(c4Key, c3Key), I6Key);
	navGraph.addEdge(I6Key, I2Key);
	//navGraph2gcs.emplace(std::make_pair(I6Key, I2Key), c3Key); //Revisar!!
	navGraph.addEdge(I2Key, I6Key);
	//navGraph2gcs.emplace(std::make_pair(I2Key, I6Key), c3Key); //Revisar!!
	navGraph.addEdge(I6Key, I3Key);
	//navGraph2gcs.emplace(std::make_pair(I6Key, I3Key), c3Key);
	navGraph.addEdge(I3Key, I6Key); 
	//navGraph2gcs.emplace(std::make_pair(I3Key, I6Key), c3Key);
	navGraph.addEdge(I6Key, I4Key);
	//navGraph2gcs.emplace(std::make_pair(I6Key, I4Key), c3Key); //Revisar!!
	navGraph.addEdge(I4Key, I6Key);
	//navGraph2gcs.emplace(std::make_pair(I5Key, I6Key), c3Key); //Revisar!!
	navGraph.addEdge(I6Key, targetKey);
	//navGraph2gcs.emplace(std::make_pair(I6Key, targetKey), c4Key);

	std::cout << "GCS Graph" << std::endl;
	gcs.print();
	std::cout << "Nav Graph" << std::endl;
	navGraph.print();
	int iters = 0;
	navGraph.printConvexSets(fout,iters);

	ConvexRelaxationMinDistanceSPP_GCS solver(&navGraph,startKey,targetKey);
	solver.optimize();
	std::cout << "Relaxed solution" << std::endl;
	std::cout << "cost: " << solver.perspectiveSolution.cost << std::endl;
	std::cout << "y: " << solver.perspectiveSolution.y << std::endl;

	std::vector<double> weights(solver.perspectiveSolution.y.size());
	Eigen::VectorXd::Map(&weights[0], solver.perspectiveSolution.y.size()) = solver.perspectiveSolution.y;
	navGraph.printGraph(fout, iters, weights);
	fout << "figure; " << std::endl;
	fout << "fGraph" << iters << "=plot(g" << iters << ",'EdgeLabel',str2num(num2str(g" << iters << ".Edges.Weights,3)));" << std::endl;
	//solver.computeFeasibleSolution(200);
	fout << "X = [];" << std::endl;
	fout << "Y = [];" << std::endl;
	fout << "for i = 1:length(centroidNavGraph" << iters << ")" << std::endl;
	fout << "	if (~isempty(centroidNavGraph" << iters << "{i}))" << std::endl;
	fout << "		X = [X; centroidNavGraph" << iters << "{i}(1)];" << std::endl;
	fout << "		Y = [Y; centroidNavGraph" << iters << "{i}(2)];" << std::endl;
	fout << "	end" << std::endl;
	fout << "end" << std::endl;
	fout << "set(fGraph" << iters << ", 'XData', X, 'YData', Y);" << std::endl;
	return 0;
}

int main(int argc, char** argv)
{
	analyse_convex_relaxation_ambiguities();
}