#include "AStarIRISConic.h"
#include "PolyhedronNode.h"
#include "PolyhedronTerminalNode.h"
#include "PointNode.h"
#include "EllipsoidNode.h"
#include "IRISConic.h"
#include "PolyhedronObstacleCircularRobotNode.h"
#include "CObsConic.h"
#include "GCS.h"
#include "ConvexRelaxationMinDistanceSolver.h"
#include "MinDistanceSolver.h"


AStarIRISConic::AStarIRISConic(CObsConic& cObs, const AStarIRISParams_t& params) : ExpandableIRISConic(cObs,params.ExpandableIRISParams), params(params)
{

}

void AStarIRISConic::do_RelaxedSolver()
{
	this->addConvexSets(qStartNode->point.p);
	bool phase1_finished = this->gcs.contains(qTargetNode->point.p);
	while (!phase1_finished)
	{
		ConvexRelaxationMinDistanceSolver solver(this->navGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
		solver.setTask();
		solver.solve();
		//solver.getGreedyPath();
		solver.computeFeasibleSolution(this->params.ExpandableIRISParams.maxItersOptimalPath);
		//Select the terminal node connected to the target that generated the feasible solution (the size of nodeKeys is always at least 2, because it includes the start and target nodes)
		int terminalNodeKey = solver.optimalPath.nodeKeys[solver.optimalPath.nodeKeys.size() - 2];
		Node* node = this->navGraph.getNode(terminalNodeKey);
		if (dynamic_cast<PolyhedronTerminalNode*>(node) == NULL)
		{
			std::cout << "We have found a feasible solution" << std::endl;
			//solver.computeFeasibleSolution();
			//this->feasibleSolution = solver.feasibleSolution;
			break;
		}
		std::cout << "Terminal node to expand " << terminalNodeKey << std::endl;
		this->expandTerminalNode(terminalNodeKey);
		/*std::cout << "Graph of convex sets" << std::endl;
		std::vector<int> nodeKeys = this->gcs.getNodeKeys();
		for (int i = 0; i < this->gcs.numNodes; i++)
		{
			std::cout << "Convex set " << nodeKeys[i] << std::endl;
			Node* node = this->gcs.getNode(nodeKeys[i]);
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			polyNode->polyhedron.print();
		}
		this->gcs.print();
		std::cout << "NavGraph of convex sets" << std::endl;
		std::vector<int> edgeNodeKeys = this->navGraph.getNodeKeys();
		for (int i = 0; i < this->navGraph.numNodes; i++)
		{
			std::cout << "Edge Convex set " << edgeNodeKeys[i] << std::endl;
			Node* node = this->navGraph.getNode(edgeNodeKeys[i]);
			if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
				PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
				polyNode->polyhedron.print();
			}
			else if (dynamic_cast<PointNode*>(node) != NULL) {
				PointNode* pointNode = (PointNode*)node->getNodeData();
				pointNode->point.print();
			}
		}
		this->navGraph.print();
		std::cout << "Relaxed Solution " << std::endl;
		std::cout << "x variables " << std::endl;
		std::cout << solver.relaxedSolution.x << std::endl;
		std::cout << "y variables " << std::endl;
		std::cout << solver.relaxedSolution.y << std::endl;
		std::cout << "z variables " << std::endl;
		std::cout << solver.relaxedSolution.z << std::endl;
		std::cout << "p variables " << std::endl;
		std::cout << solver.relaxedSolution.p << std::endl;
		std::cout << "l variables " << std::endl;
		std::cout << solver.relaxedSolution.l << std::endl;
		std::cout << "cost " << std::endl;
		std::cout << solver.relaxedSolution.cost << std::endl;*/
		this->optimalPath = solver.optimalPath;
		/*std::cout << "Feasible Solution " << std::endl;
		std::cout << "x variables " << std::endl;
		std::cout << solver.feasibleSolution.x << std::endl;
		std::cout << "nodes keys ";
		for (std::vector<int>::iterator it = solver.feasibleSolution.nodeKeys.begin(); it != solver.feasibleSolution.nodeKeys.end(); it++)
			std::cout << *it << " ";
		std::cout << std::endl;
		std::cout << "edge keys ";
		for (std::vector<int>::iterator it = solver.feasibleSolution.edgeKeys.begin(); it != solver.feasibleSolution.edgeKeys.end(); it++)
			std::cout << *it << " ";
		std::cout << std::endl;
		std::cout << "cost " << std::endl;
		std::cout << solver.feasibleSolution.cost << std::endl;*/
		/*phase1_finished = true;
		std::vector<int> nodeKeys = this->navGraph.getNodeKeys();
		for (std::vector<int>::iterator it = nodeKeys.begin(); it != nodeKeys.end(); it++)
		{
			Node* node = (Node*)this->navGraph.getNode(*it);
			if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
			{
				phase1_finished = false;
				break;
			}
		}*/
		phase1_finished = this->gcs.contains(qTargetNode->point.p);
		//if (phase1_finished)
		//{
		//	solver.computeFeasibleSolution();
		//	this->feasibleSolution = solver.feasibleSolution;
		//}
	}
}



void AStarIRISConic::do_RelaxedSolver(std::ostream& out)
{
	bool dumpResults = true;
	int iters = 1;
	std::cout << "Starting Phase 1" << std::endl;
	if (dumpResults)
	{
		out << "% Starting Phase 1" << std::endl;
		out << "Frames=[];" << std::endl;
	}
	this->addConvexSet(qStartNode->point.p);
	bool phase1_finished = this->gcs.contains(qTargetNode->point.p);
	while (!phase1_finished)
	{
		ConvexRelaxationMinDistanceSolver solver(this->navGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
		solver.setTask();
		solver.solve();
		//Path_t nodePath = solver.getGreedyPath();
		//Path_t nodePath = solver.getMCPath();  //Just one iteration
		solver.computeFeasibleSolution(this->params.ExpandableIRISParams.maxItersOptimalPath);
		Path_t nodePath = solver.optimalPath;

		if (dumpResults)
		{
			std::vector<double> weights(solver.relaxedSolution.y.size());
			Eigen::VectorXd::Map(&weights[0], solver.relaxedSolution.y.size()) = solver.relaxedSolution.y;
			this->navGraph.printConvexSets(out, iters); 
			this->navGraph.printGraph(out, iters, weights);

			out << "figure; " << std::endl;
			//out << "fGraph" << iters << "=plot(g" << iters << ",'EdgeLabel',str2num(num2str(g" << iters << ".Edges.Weights,3)));" << std::endl;
			out << "fGraph" << iters << "=plot(g" << iters << ");" << std::endl;


			out << "nodesGraphPath" << iters << "={";
			for (std::vector<int>::iterator it = nodePath.nodeKeys.begin(); it != nodePath.nodeKeys.end(); it++)
			{
				out << "'" << (*it) + 1 << "' ";
			}
			out << "}';" << std::endl;

			out << "highlight(fGraph" << iters << ",nodesGraphPath" << iters << ",'NodeColor','r','EdgeColor','r','LineWidth',3)" << std::endl;
		}
		int terminalNodeKey = nodePath.nodeKeys[nodePath.nodeKeys.size() - 2];
		Node* node = this->navGraph.getNode(terminalNodeKey);
		if (dynamic_cast<PolyhedronTerminalNode*>(node) == NULL)
		{
			std::cout << "We have found a feasible solution using intersection node key " << terminalNodeKey << std::endl;
			//solver.computeFeasibleSolution();
			//this->feasibleSolution = solver.feasibleSolution;
			break;
		}
		
		if (dumpResults)
		{
			this->navGraph.printConvexSet(out, iters,terminalNodeKey);
			this->gcs.print(out, iters);

			out << "allA=[AObs AGCS" << iters << " ANavGraph" << iters << "]; " << std::endl;
			out << "allB=[bObs bGCS" << iters << " bNavGraph" << iters << "];" << std::endl;
			out << "allColors=[colObs colGCS" << iters << " colNavGraph" << iters << "];" << std::endl;
			out << "fGCS" << iters << "=figure; " << std::endl;
			out << "plotregion(allA, allB, [], [], allColors);" << std::endl;
			out << "hold on;" << std::endl;
			out << "xlabel('X [m.]');" << std::endl;
			out << "ylabel('Y [m.]');" << std::endl;
			out << "for i = 1:length(centroidGCS" << iters << ")" << std::endl;
			out << "  text(centroidGCS" << iters << "{i}(1), centroidGCS" << iters << "{i}(2),nameGCS" << iters << "{i}, 'FontSize', 12); " << std::endl;
			out << "end" << std::endl;
			out << "for i = 1:length(centroidNavGraph" << iters << ")" << std::endl;
			out << "  if (~isempty(centroidNavGraph" << iters << "{i}))" << std::endl;
			out << "    text(centroidNavGraph" << iters << "{i}(1), centroidNavGraph" << iters << "{i}(2),nameNavGraph" << iters << "{i},'FontSize',6,'Color','k');" << std::endl;
			out << "  end" << std::endl; 
			out << "end" << std::endl;
			out << "axis equal;" << std::endl;
			out << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;
			out << "Frames=[Frames;getframe(gcf)];" << std::endl;
		}

		std::cout << "Terminal node to expand " << terminalNodeKey << std::endl;
		this->expandTerminalNode(terminalNodeKey);
		phase1_finished = this->gcs.contains(qTargetNode->point.p);
		iters++;
	}

	std::cout << "Phase 1 finished!" << std::endl;

	std::cout << "Starting Phase 2" << std::endl;
	if (dumpResults)
		out << "% Starting Phase 2" << std::endl;
	//Now compute a simplified graph to obtain a feasible solution
	NavGraph simplifiedGraphBest = this->getGraphWithoutTerminalConnections();

	ConvexRelaxationMinDistanceSolver solverTargetReached(simplifiedGraphBest, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
	solverTargetReached.setTask();
	solverTargetReached.solve();
	solverTargetReached.computeFeasibleSolution(this->params.ExpandableIRISParams.maxItersOptimalPath);
	double bestCost = solverTargetReached.feasibleSolution.cost;
	Eigen::MatrixXd optimalTargetReachedPath= solverTargetReached.feasibleSolution.x;
	Path_t nodesOptimalTargetReachedGraphPath = solverTargetReached.optimalPath;

	std::cout << "Initial best cost " << bestCost << std::endl;
	
	if (dumpResults)
	{
		{
			this->navGraph.printConvexSets(out, iters);
			simplifiedGraphBest.printGraph(out, iters);
			//this->navGraph.printGraph(out, iters);

			out << "optimalPath" << iters << "=[";
			for (int i = 0; i < solverTargetReached.feasibleSolution.x.rows(); i++)
			{
				for (int j = 0; j < solverTargetReached.feasibleSolution.x.cols(); j++)
				{
					if (j < (solverTargetReached.feasibleSolution.x.cols() - 1))
						out << solverTargetReached.feasibleSolution.x(i, j) << " ";
					else
						out << solverTargetReached.feasibleSolution.x(i, j);
				}
				if (i < (solverTargetReached.feasibleSolution.x.rows() - 1))
					out << ";";
			}
			out << "];" << std::endl;
			

			out << "nodesOptimalGraphPath" << iters << "={";
			for (std::vector<int>::iterator it = nodesOptimalTargetReachedGraphPath.nodeKeys.begin(); it != nodesOptimalTargetReachedGraphPath.nodeKeys.end(); it++)
			{
				out << "'" << (*it) + 1 << "' ";
			}
			out << "}';" << std::endl;
			
			

			out << "figure; " << std::endl;
			out << "fGraph" << iters << "=plot(g" << iters << ");" << std::endl;


			out << "nodesGraphPath" << iters << "={";
			for (std::vector<int>::iterator it = nodesOptimalTargetReachedGraphPath.nodeKeys.begin(); it != nodesOptimalTargetReachedGraphPath.nodeKeys.end(); it++)
			{
				out << "'" << (*it) + 1 << "' ";
			}
			out << "}';" << std::endl;

			out << "highlight(fGraph" << iters << ",nodesGraphPath" << iters << ",'NodeColor','r','EdgeColor','r','LineWidth',3)" << std::endl;
			
			this->gcs.print(out, iters);

			out << "allA=[AObs AGCS" << iters << "]; " << std::endl;
			out << "allB=[bObs bGCS" << iters << "];" << std::endl;
			out << "allColors=[colObs colGCS" << iters << "];" << std::endl;
			out << "allColors=[colObs colGCS" << iters << "];" << std::endl;
			out << "fGCS" << iters << "=figure; " << std::endl;
			out << "plotregion(allA, allB, [], [], allColors);" << std::endl;
			out << "hold on;" << std::endl;
			out << "xlabel('X [m.]');" << std::endl;
			out << "ylabel('Y [m.]');" << std::endl;
			out << "for i = 1:length(centroidGCS" << iters << ")" << std::endl;
			out << "  text(centroidGCS" << iters << "{i}(1), centroidGCS" << iters << "{i}(2),nameGCS" << iters << "{i}, 'FontSize', 12); " << std::endl;
			out << "end" << std::endl;
			out << "for i = 1:length(centroidNavGraph" << iters << ")" << std::endl;
			out << "  if (~isempty(centroidNavGraph" << iters << "{i}))" << std::endl;
			out << "    text(centroidNavGraph" << iters << "{i}(1), centroidNavGraph" << iters << "{i}(2),nameNavGraph" << iters << "{i},'FontSize',6,'Color','k');" << std::endl;
			out << "  end" << std::endl;
			out << "end" << std::endl;
			out << "axis equal;" << std::endl;
			out << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;
			out << "plot(optimalPath" << iters << "(:,1),optimalPath" << iters << "(:,2),'b','LineWidth',3);";
			out << "Frames=[Frames;getframe(gcf)];" << std::endl;
		}
		iters++;
	}

	//return;

	int count = 0;
	while (true)
	{
		count++;
		std::cout << "Count: " << count << std::endl;
		//if (count > 3)
		//	break;
		//Obtain the optimal and relaxed solutions againg of the current navGraph to plot the results
		ConvexRelaxationMinDistanceSolver solverPhase2(this->navGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
		solverPhase2.setTask();
		solverPhase2.solve();
		solverPhase2.computeFeasibleSolution(this->params.ExpandableIRISParams.maxItersOptimalPath);


		if (solverPhase2.feasibleSolution.cost >= bestCost)
			break;

		int terminalNodeKey = solverPhase2.optimalPath.nodeKeys[solverPhase2.optimalPath.nodeKeys.size()-2];
		//if (count == 1)
		//	terminalNodeKey = 21;
		Node* node = this->navGraph.getNode(terminalNodeKey);

		if (dumpResults)
		{
			this->navGraph.printConvexSet(out, iters, terminalNodeKey); //We need to do this before expanding, because later the node will be removed
		}

		if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
		{
		}
		else
		{
			optimalTargetReachedPath = solverPhase2.feasibleSolution.x;
			nodesOptimalTargetReachedGraphPath = solverPhase2.optimalPath;
			bestCost = solverPhase2.feasibleSolution.cost;
			std::cout << "Optimal cost " << bestCost << std::endl;
			break;
		}

		if (dumpResults)
		{
			this->navGraph.printConvexSets(out, iters);
			this->navGraph.printGraph(out, iters);
			out << "optimalPath" << iters << "=[";
			for (int i = 0; i < optimalTargetReachedPath.rows(); i++)
			{
				for (int j = 0; j < optimalTargetReachedPath.cols(); j++)
				{
					if (j < (optimalTargetReachedPath.cols() - 1))
						out << optimalTargetReachedPath(i, j) << " ";
					else
						out << optimalTargetReachedPath(i, j);
				}
				if (i < (optimalTargetReachedPath.rows() - 1))
					out << ";";
			}
			out << "];" << std::endl;

			if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
			{
				out << "optimisticPath" << iters << "=[";
				for (int i = 0; i < solverPhase2.feasibleSolution.x.rows(); i++)
				{
					for (int j = 0; j < solverPhase2.feasibleSolution.x.cols(); j++)
					{
						if (j < (solverPhase2.feasibleSolution.x.cols() - 1))
							out << solverPhase2.feasibleSolution.x(i, j) << " ";
						else
							out << solverPhase2.feasibleSolution.x(i, j);
					}
					if (i < (solverPhase2.feasibleSolution.x.rows() - 1))
						out << ";";
				}
				out << "];" << std::endl;
			}

			out << "nodesOptimalGraphPath" << iters << "={";
			for (std::vector<int>::iterator it = nodesOptimalTargetReachedGraphPath.nodeKeys.begin(); it != nodesOptimalTargetReachedGraphPath.nodeKeys.end(); it++)
			{
				out << "'" << (*it) + 1 << "' ";
			}
			out << "}';" << std::endl;

			out << "figure; " << std::endl;
			out << "fGraph" << iters << "=plot(g" << iters << ");" << std::endl;

			out << "highlight(fGraph" << iters << ",nodesOptimalGraphPath" << iters << ",'NodeColor','r','EdgeColor','g','LineWidth',1.5)" << std::endl;


			this->gcs.print(out, iters);

			out << "allA=[AObs AGCS" << iters << " ANavGraph" << iters << "]; " << std::endl;
			out << "allB=[bObs bGCS" << iters << " bNavGraph" << iters << "];" << std::endl;
			out << "allColors=[colObs colGCS" << iters << " colNavGraph" << iters << "];" << std::endl;
			out << "fGCS" << iters << "=figure; " << std::endl;
			out << "plotregion(allA, allB, [], [], allColors);" << std::endl;
			out << "hold on;" << std::endl;
			out << "xlabel('X [m.]');" << std::endl;
			out << "ylabel('Y [m.]');" << std::endl;
			out << "for i = 1:length(centroidGCS" << iters << ")" << std::endl;
			out << "  text(centroidGCS" << iters << "{i}(1), centroidGCS" << iters << "{i}(2),nameGCS" << iters << "{i},'FontSize',12);" << std::endl;
			out << "end" << std::endl;
			out << "for i = 1:length(centroidNavGraph" << iters << ")" << std::endl;
			out << "  if (~isempty(centroidNavGraph" << iters << "{i}))" << std::endl;
			out << "    text(centroidNavGraph" << iters << "{i}(1), centroidNavGraph" << iters << "{i}(2),nameNavGraph" << iters << "{i},'FontSize',6,'Color','k');" << std::endl;
			out << "  end" << std::endl;
			out << "end" << std::endl;
			out << "axis equal;" << std::endl;
			out << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;
			if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
			{
				out << "plot(optimisticPath" << iters << "(:,1),optimisticPath" << iters << "(:,2),'b.','LineWidth',3);";
			}
			out << "plot(optimalPath" << iters << "(:,1),optimalPath" << iters << "(:,2),'b','LineWidth',3);";
			out << "Frames=[Frames;getframe(gcf)];" << std::endl;
			
		}
		iters++;

		if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
		{
			std::cout << "Terminal node to expand " << terminalNodeKey << std::endl;
			this->expandTerminalNode(terminalNodeKey);
		}
	}

	std::cout << "Phase 2 finished!" << std::endl;
	std::cout << "Now it's computing the optimal path with the final GCS graph (just for dumping the results or it uses the MIP solver)" << std::endl;
	//Compute one more time to obtain the optimal path with the resulting graph
	simplifiedGraphBest = this->getGraphWithoutTerminalConnections();
	ConvexRelaxationMinDistanceSolver solverFinished(simplifiedGraphBest, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
	//MinDistanceSolver solverFinished(simplifiedGraphBest, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
	solverFinished.setTask();
	solverFinished.solve();
	solverFinished.computeFeasibleSolution(10*this->params.ExpandableIRISParams.maxItersOptimalPath);
	bestCost = solverFinished.feasibleSolution.cost;
	optimalTargetReachedPath = solverFinished.feasibleSolution.x;
	nodesOptimalTargetReachedGraphPath = solverFinished.optimalPath;
	std::cout << "Optimal cost " << bestCost << std::endl;
	std::cout << "Finished!" << std::endl;

	if (dumpResults)
	{
		{
			simplifiedGraphBest.printConvexSets(out, iters);
			simplifiedGraphBest.printGraph(out, iters);
			//this->navGraph.printGraph(out, iters);


			out << "optimalPath" << iters << "=[";
			for (int i = 0; i < optimalTargetReachedPath.rows(); i++)
			{
				for (int j = 0; j < optimalTargetReachedPath.cols(); j++)
				{
					if (j < (optimalTargetReachedPath.cols() - 1))
						out << optimalTargetReachedPath(i, j) << " ";
					else
						out << optimalTargetReachedPath(i, j);
				}
				if (i < (optimalTargetReachedPath.rows() - 1))
					out << ";";
			}
			out << "];" << std::endl;
			out << "nodesOptimalGraphPath" << iters << "={";
			for (std::vector<int>::iterator it = nodesOptimalTargetReachedGraphPath.nodeKeys.begin(); it != nodesOptimalTargetReachedGraphPath.nodeKeys.end(); it++)
			{
				out << "'" << (*it) + 1 << "' ";
			}
			out << "}';" << std::endl;



			out << "figure; " << std::endl;
			out << "fGraph" << iters << "=plot(g" << iters << ");" << std::endl;


			out << "nodesGraphPath" << iters << "={";
			for (std::vector<int>::iterator it = nodesOptimalTargetReachedGraphPath.nodeKeys.begin(); it != nodesOptimalTargetReachedGraphPath.nodeKeys.end(); it++)
			{
				out << "'" << (*it) + 1 << "' ";
			}
			out << "}';" << std::endl;

			out << "highlight(fGraph" << iters << ",nodesGraphPath" << iters << ",'NodeColor','r','EdgeColor','r','LineWidth',3)" << std::endl;

			this->gcs.print(out, iters);

			out << "allA=[AObs AGCS" << iters << "]; " << std::endl;
			out << "allB=[bObs bGCS" << iters << "];" << std::endl;
			out << "allColors=[colObs colGCS" << iters << "];" << std::endl;
			out << "allColors=[colObs colGCS" << iters << "];" << std::endl;
			out << "fGCS" << iters << "=figure; " << std::endl;
			out << "plotregion(allA, allB, [], [], allColors);" << std::endl;
			out << "hold on;" << std::endl;
			out << "xlabel('X [m.]');" << std::endl;
			out << "ylabel('Y [m.]');" << std::endl;
			out << "for i = 1:length(centroidGCS" << iters << ")" << std::endl;
			out << "  text(centroidGCS" << iters << "{i}(1), centroidGCS" << iters << "{i}(2),nameGCS" << iters << "{i}, 'FontSize', 12); " << std::endl;
			out << "end" << std::endl;
			out << "for i = 1:length(centroidNavGraph" << iters << ")" << std::endl;
			out << "  if (~isempty(centroidNavGraph" << iters << "{i}))" << std::endl;
			out << "    text(centroidNavGraph" << iters << "{i}(1), centroidNavGraph" << iters << "{i}(2),nameNavGraph" << iters << "{i},'FontSize',6,'Color','k');" << std::endl;
			out << "  end" << std::endl;
			out << "end" << std::endl;
			out << "axis equal;" << std::endl;
			out << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;
			out << "plot(optimalPath" << iters << "(:,1),optimalPath" << iters << "(:,2),'b','LineWidth',3);";
			out << "Frames=[Frames;getframe(gcf)];" << std::endl;
		}
		iters++;
	}
}

void AStarIRISConic::do_MIPSolver(std::ostream& out)
{
	{
		bool dumpResults = true;
		int iters = 1;
		std::cout << "Starting Phase 1" << std::endl;
		this->addConvexSets(qStartNode->point.p);
		bool phase1_finished = this->gcs.contains(qTargetNode->point.p);
		while (!phase1_finished)
		{
			MinDistanceSolver solverMIP(this->navGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
			solverMIP.setTask();
			solverMIP.solve();
			solverMIP.computeFeasibleSolution();
			Path_t optimalPath = solverMIP.optimalPath;

			if (dumpResults)
			{
				this->navGraph.printConvexSets(out, iters);
				this->navGraph.printGraph(out, iters);
				out << "optimalPath" << iters << "=[";
				for (int i = 0; i < solverMIP.feasibleSolution.x.rows(); i++)
				{
					for (int j = 0; j < solverMIP.feasibleSolution.x.cols(); j++)
					{
						if (j < (solverMIP.feasibleSolution.x.cols() - 1))
							out << solverMIP.feasibleSolution.x(i, j) << " ";
						else
							out << solverMIP.feasibleSolution.x(i, j);
					}
					if (i < (solverMIP.feasibleSolution.x.rows() - 1))
						out << ";";
				}
				out << "];" << std::endl;

				out << "nodesOptimalGraphPath" << iters << "={";
				for (std::vector<int>::iterator it = optimalPath.nodeKeys.begin(); it != optimalPath.nodeKeys.end(); it++)
				{
					out << "'" << (*it) + 1 << "' ";
				}
				out << "}';" << std::endl;

				out << "figure; " << std::endl;
				out << "fGraph" << iters << "=plot(g" << iters << ");" << std::endl;

				out << "highlight(fGraph" << iters << ",nodesOptimalGraphPath" << iters << ",'NodeColor','r','EdgeColor','g','LineWidth',1.5)" << std::endl;
				//Select the terminal node connected to the target that generated the feasible solution (the size of nodeKeys is always at least 2, because it includes the start and target nodes)
			}
			
			int terminalNodeKey = optimalPath.nodeKeys[optimalPath.nodeKeys.size() - 2];
			

			Node* node = this->navGraph.getNode(terminalNodeKey);
			if (dynamic_cast<PolyhedronTerminalNode*>(node) == NULL)
			{
				std::cout << "We have found a feasible solution" << std::endl;
				//solver.computeFeasibleSolution();
				//this->feasibleSolution = solver.feasibleSolution;
				break;
			}
			if (dumpResults)
			{
				this->navGraph.printConvexSet(out, iters, terminalNodeKey);
				this->gcs.print(out, iters);
				
				out << "allA=[AObs AGCS" << iters << " ANavGraph" << iters << "]; " << std::endl;
				out << "allB=[bObs bGCS" << iters << " bNavGraph" << iters << "];" << std::endl;
				out << "allColors=[colObs colGCS" << iters << " colNavGraph" << iters << "];" << std::endl;
				out << "fGCS" << iters << "=figure; " << std::endl;
				out << "plotregion(allA, allB, [], [], allColors);" << std::endl;
				out << "hold on;" << std::endl;
				out << "xlabel('X [m.]');" << std::endl;
				out << "ylabel('Y [m.]');" << std::endl;
				out << "for i = 1:length(centroidGCS" << iters << ")" << std::endl;
				out << "  text(centroidGCS" << iters << "{i}(1), centroidGCS" << iters << "{i}(2),nameGCS" << iters << "{i}, 'FontSize', 12); " << std::endl;
				out << "end" << std::endl;
				out << "for i = 1:length(centroidNavGraph" << iters << ")" << std::endl;
				out << "  text(centroidNavGraph" << iters << "{i}(1), centroidNavGraph" << iters << "{i}(2),nameNavGraph" << iters << "{i},'FontSize',6,'Color','k');" << std::endl;
				out << "end" << std::endl;
				out << "axis equal;" << std::endl;
				out << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;

				//out << "plot(centroidNavGraph" << iters << "{" << (idx+1) << "}(1), centroidNavGraph" << iters << "{" << (idx+1) << "}(2),'ro','MarkerSize',12,'LineWidth',3);" << std::endl;
				out << "plot(optimalPath" << iters << "(:,1),optimalPath" << iters << "(:,2),'b.','LineWidth',3);";

			}

			std::cout << "Terminal node to expand " << terminalNodeKey << std::endl;
			
			this->expandTerminalNode(terminalNodeKey);
			phase1_finished = this->gcs.contains(qTargetNode->point.p);
			iters++;
		}


		std::cout << "Starting Phase 2" << std::endl;
		//Now compute a simplified graph to obtain a feasible solution
		NavGraph simplifiedGraphBest = this->getGraphWithoutTerminalConnections();
		MinDistanceSolver solverMIPTargetReached(simplifiedGraphBest, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
		solverMIPTargetReached.setTask();
		solverMIPTargetReached.solve();
		solverMIPTargetReached.computeFeasibleSolution();
		double bestCost = solverMIPTargetReached.feasibleSolution.cost;
		Eigen::MatrixXd optimalTargetReachedPath = solverMIPTargetReached.feasibleSolution.x;
		Path_t nodesOptimalTargetReachedGraphPath = solverMIPTargetReached.optimalPath;

		if (dumpResults)
		{
			this->navGraph.printConvexSets(out, iters);
			simplifiedGraphBest.printGraph(out, iters);
			out << "optimalPath" << iters << "=[";
			for (int i = 0; i < solverMIPTargetReached.feasibleSolution.x.rows(); i++)
			{
				for (int j = 0; j < solverMIPTargetReached.feasibleSolution.x.cols(); j++)
				{
					if (j < (solverMIPTargetReached.feasibleSolution.x.cols() - 1))
						out << solverMIPTargetReached.feasibleSolution.x(i, j) << " ";
					else
						out << solverMIPTargetReached.feasibleSolution.x(i, j);
				}
				if (i < (solverMIPTargetReached.feasibleSolution.x.rows() - 1))
					out << ";";
			}
			out << "];" << std::endl;

			out << "nodesOptimalGraphPath" << iters << "={";
			for (std::vector<int>::iterator it = nodesOptimalTargetReachedGraphPath.nodeKeys.begin(); it != nodesOptimalTargetReachedGraphPath.nodeKeys.end(); it++)
			{
				out << "'" << (*it) + 1 << "' ";
			}
			out << "}';" << std::endl;

			out << "figure; " << std::endl;
			out << "fGraph" << iters << "=plot(g" << iters << ");" << std::endl;

			out << "highlight(fGraph" << iters << ",nodesOptimalGraphPath" << iters << ",'NodeColor','r','EdgeColor','g','LineWidth',1.5)" << std::endl;

			this->gcs.print(out, iters);
			
			out << "allA=[AObs AGCS" << iters << "]; " << std::endl;
			out << "allB=[bObs bGCS" << iters << "];" << std::endl;
			out << "allColors=[colObs colGCS" << iters << "];" << std::endl;
			out << "fGCS" << iters << "=figure; " << std::endl;
			out << "plotregion(allA, allB, [], [], allColors);" << std::endl;
			out << "hold on;" << std::endl;
			out << "xlabel('X [m.]');" << std::endl;
			out << "ylabel('Y [m.]');" << std::endl;
			out << "for i = 1:length(centroidGCS" << iters << ")" << std::endl;
			out << "  text(centroidGCS" << iters << "{i}(1), centroidGCS" << iters << "{i}(2),nameGCS" << iters << "{i}, 'FontSize', 12); " << std::endl;
			out << "end" << std::endl;
			out << "for i = 1:length(centroidNavGraph" << iters << ")" << std::endl;
			out << "  text(centroidNavGraph" << iters << "{i}(1), centroidNavGraph" << iters << "{i}(2),nameNavGraph" << iters << "{i},'FontSize',6,'Color','k');" << std::endl;
			out << "end" << std::endl;
			out << "axis equal;" << std::endl;
			out << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;
			out << "plot(optimalPath" << iters << "(:,1),optimalPath" << iters << "(:,2),'b','LineWidth',3);";
			
		}
		iters++;

		int count = 0;
		while (true)
		{
			count++;
			std::cout << "Count: " << count << std::endl;
			//if (count > 3)
			//	break;
			//Obtain the optimal and relaxed solutions againg of the current navGraph to plot the results
			MinDistanceSolver solverMIPPhase2(this->navGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
			solverMIPPhase2.setTask();
			solverMIPPhase2.solve();
			solverMIPPhase2.computeFeasibleSolution();
			
			if (solverMIPPhase2.feasibleSolution.cost > bestCost)
				break;

			int terminalNodeKey = solverMIPPhase2.optimalPath.nodeKeys[solverMIPPhase2.optimalPath.nodeKeys.size()-2];
			Node* node = this->navGraph.getNode(terminalNodeKey);
			
			if (dumpResults)
			{
				this->navGraph.printConvexSet(out, iters, terminalNodeKey); //We need to do this before expanding, because later the node will be removed
			}

			if (dynamic_cast<PolyhedronNode*>(node) != NULL)
			{
				optimalTargetReachedPath = solverMIPPhase2.feasibleSolution.x;
				nodesOptimalTargetReachedGraphPath = solverMIPPhase2.optimalPath;
				bestCost = solverMIPPhase2.feasibleSolution.cost;
			}

			if (dumpResults)
			{
				this->navGraph.printConvexSets(out, iters);
				this->navGraph.printGraph(out, iters);
				out << "optimalPath" << iters << "=[";
				for (int i = 0; i < optimalTargetReachedPath.rows(); i++)
				{
					for (int j = 0; j < optimalTargetReachedPath.cols(); j++)
					{
						if (j < (optimalTargetReachedPath.cols() - 1))
							out << optimalTargetReachedPath(i, j) << " ";
						else
							out << optimalTargetReachedPath(i, j);
					}
					if (i < (optimalTargetReachedPath.rows() - 1))
						out << ";";
				}
				out << "];" << std::endl;

				if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
				{
					out << "optimisticPath" << iters << "=[";
					for (int i = 0; i < solverMIPPhase2.feasibleSolution.x.rows(); i++)
					{
						for (int j = 0; j < solverMIPPhase2.feasibleSolution.x.cols(); j++)
						{
							if (j < (solverMIPPhase2.feasibleSolution.x.cols() - 1))
								out << solverMIPPhase2.feasibleSolution.x(i, j) << " ";
							else
								out << solverMIPPhase2.feasibleSolution.x(i, j);
						}
						if (i < (solverMIPPhase2.feasibleSolution.x.rows() - 1))
							out << ";";
					}
					out << "];" << std::endl;
				}

				out << "nodesOptimalGraphPath" << iters << "={";
				for (std::vector<int>::iterator it = nodesOptimalTargetReachedGraphPath.nodeKeys.begin(); it != nodesOptimalTargetReachedGraphPath.nodeKeys.end(); it++)
				{
					out << "'" << (*it) + 1 << "' ";
				}
				out << "}';" << std::endl;

				out << "figure; " << std::endl;
				out << "fGraph" << iters << "=plot(g" << iters << ");" << std::endl;

				out << "highlight(fGraph" << iters << ",nodesOptimalGraphPath" << iters << ",'NodeColor','r','EdgeColor','g','LineWidth',1.5)" << std::endl;


				this->gcs.print(out, iters);

				out << "allA=[AObs AGCS" << iters << " ANavGraph" << iters << "]; " << std::endl;
				out << "allB=[bObs bGCS" << iters << " bNavGraph" << iters << "];" << std::endl;
				out << "allColors=[colObs colGCS" << iters << " colNavGraph" << iters << "];" << std::endl;
				out << "fGCS" << iters << "=figure; " << std::endl;
				out << "plotregion(allA, allB, [], [], allColors);" << std::endl;
				out << "hold on;" << std::endl;
				out << "xlabel('X [m.]');" << std::endl;
				out << "ylabel('Y [m.]');" << std::endl;
				out << "for i = 1:length(centroidGCS" << iters << ")" << std::endl;
				out << "  text(centroidGCS" << iters << "{i}(1), centroidGCS" << iters << "{i}(2),nameGCS" << iters << "{i},'FontSize',12);" << std::endl;
				out << "end" << std::endl;
				out << "for i = 1:length(centroidNavGraph" << iters << ")" << std::endl;
				out << "  text(centroidNavGraph" << iters << "{i}(1), centroidNavGraph" << iters << "{i}(2),nameNavGraph" << iters << "{i},'FontSize',6,'Color','k');" << std::endl;
				out << "end" << std::endl;
				out << "axis equal;" << std::endl;
				out << "axis(reshape(Range.range,1,numel(Range.range)));" << std::endl;
				if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
				{
					out << "plot(optimisticPath" << iters << "(:,1),optimisticPath" << iters << "(:,2),'b.','LineWidth',3);";
				}
				out << "plot(optimalPath" << iters << "(:,1),optimalPath" << iters << "(:,2),'b','LineWidth',3);";
				iters++;
			}
			
			if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
			{
				std::cout << "Terminal node to expand " << terminalNodeKey << std::endl;	
				if (terminalNodeKey == 54)
					std::cout << "This one is problematic" << std::endl;
				this->expandTerminalNode(terminalNodeKey);	
			}
		}
	}
}

void AStarIRISConic::do_MIPSolver()
{
	this->do_MIPSolver(std::cout);
	this->addConvexSets(qStartNode->point.p);
	bool phase1_finished = this->gcs.contains(qTargetNode->point.p);
	while (!phase1_finished)
	{
		MinDistanceSolver solver(this->navGraph, this->qStartNodeNavGraphKey, this->qTargetNodeNavGraphKey);
		solver.setTask();
		solver.solve();
		//Select the terminal node connected to the target that generated the feasible solution (the size of nodeKeys is always at least 2, because it includes the start and target nodes)
		int terminalNodeKey = solver.optimalPath.nodeKeys[solver.optimalPath.nodeKeys.size() - 2];
		std::cout << "Terminal node to expand " << terminalNodeKey << std::endl;
		//if (terminalNodeKey == 38)
		//{
		//	std::cout << "This is the one" << std::endl; //For some strange reason the seed generator is not working properly, even if I increase the number of maximum trials, but apparently the facet is expandable
		//	//params.maxTrialsTerminal = 999;
		//}
		this->expandTerminalNode(terminalNodeKey);
		

		/*std::vector<int> nodeKeysFound = this->gcs.findConvexSets(qTargetNode->point.p);
		bool qTargetInGCS = nodeKeysFound.size()>0;
		if (qTargetInGCS)
		{
			std::cout << "Graph of convex sets" << std::endl;
			std::vector<int> nodeKeys = this->gcs.getNodeKeys();
			for (int i = 0; i < this->gcs.numNodes; i++)
			{
				std::cout << "Convex set " << nodeKeys[i] << std::endl;
				Node* node = this->gcs.getNode(nodeKeys[i]);
				PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
				polyNode->polyhedron.print();
			}
			this->gcs.print();
			std::cout << "NavGraph of convex sets" << std::endl;
			std::vector<int> edgeNodeKeys = this->navGraph.getNodeKeys();
			for (int i = 0; i < this->navGraph.numNodes; i++)
			{
				Node* node = this->navGraph.getNode(edgeNodeKeys[i]);
				if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL) {
					std::cout << "Terminal Node Convex set " << edgeNodeKeys[i] << std::endl;
					PolyhedronTerminalNode* polyNode = (PolyhedronTerminalNode*)node->getNodeData();
					polyNode->polyhedron.print();
				}
				else if (dynamic_cast<PolyhedronNode*>(node) != NULL) {
					std::cout << "Polyhedron Node Convex set " << edgeNodeKeys[i] << std::endl;
					PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
					polyNode->polyhedron.print();
				}
				else if (dynamic_cast<PointNode*>(node) != NULL) {
					std::cout << "Point Node Convex set " << edgeNodeKeys[i] << std::endl;
					PointNode* pointNode = (PointNode*)node->getNodeData();
					pointNode->point.print();
				}
			}
			this->navGraph.print();
			std::cout << "MIP Solution " << std::endl;
			std::cout << "x variables " << std::endl;
			std::cout << solver.MIPSolution.x << std::endl;
			std::cout << "y variables " << std::endl;
			std::cout << solver.MIPSolution.y << std::endl;
			std::cout << "z variables " << std::endl;
			std::cout << solver.MIPSolution.z << std::endl;
			std::cout << "p variables " << std::endl;
			std::cout << solver.MIPSolution.p << std::endl;
			std::cout << "l variables " << std::endl;
			std::cout << solver.MIPSolution.l << std::endl;
			std::cout << "cost " << std::endl;
			std::cout << solver.MIPSolution.cost << std::endl;
			this->feasibleSolution = solver.feasibleSolution;
			std::cout << "Feasible Solution " << std::endl;
			std::cout << "x variables " << std::endl;
			std::cout << solver.feasibleSolution.x << std::endl;
			std::cout << "nodes keys ";
			for (std::vector<int>::iterator it = solver.feasibleSolution.nodeKeys.begin(); it != solver.feasibleSolution.nodeKeys.end(); it++)
				std::cout << *it << " ";
			std::cout << std::endl;
			std::cout << "edge keys ";
			for (std::vector<int>::iterator it = solver.feasibleSolution.edgeKeys.begin(); it != solver.feasibleSolution.edgeKeys.end(); it++)
				std::cout << *it << " ";
			std::cout << std::endl;
			std::cout << "cost " << std::endl;
			std::cout << solver.feasibleSolution.cost << std::endl;
		}
		Eigen::VectorXd previousConf=solver.feasibleSolution.x.row(solver.feasibleSolution.x.rows()-2);
		bool previousConfValid = false;
		for (std::vector<int>::iterator itNodeKeysFound = nodeKeysFound.begin(); itNodeKeysFound != nodeKeysFound.end(); itNodeKeysFound++)
		{
			Node* node=(Node*)this->gcs.getNode(*itNodeKeysFound);
			PolyhedronNode* polyNode = (PolyhedronNode*)node->getNodeData();
			if (polyNode->polyhedron.isInside(previousConf))
			{
				previousConfValid = true;
				break;
			}
		}
		phase1_finished = qTargetInGCS && previousConfValid;*/
		/*phase1_finished = true;
		std::vector<int> nodeKeys = this->navGraph.getNodeKeys();
		for (std::vector<int>::iterator it = nodeKeys.begin(); it != nodeKeys.end(); it++)
		{
			Node* node = (Node*)this->navGraph.getNode(*it);
			if (dynamic_cast<PolyhedronTerminalNode*>(node) != NULL)
			{
				phase1_finished = false;
				break;
			}
		}*/
		phase1_finished = this->gcs.contains(qTargetNode->point.p);
	}
}


AStarIRISParams_t AStarIRISConic::getDefaultAStarIRISParams()
{
	AStarIRISParams_t params = { ExpandableIRISConic::getDefaultExpandableIRISParams()};
	return params;
}