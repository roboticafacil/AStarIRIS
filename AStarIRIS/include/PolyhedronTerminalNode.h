#ifndef POLYHEDRON_TERMINAL_NODE
#define POLYHEDRON_TERMINAL_NODE
#include<Eigen/Dense>
#include "PolyhedronNode.h"
#include "Node.h"

class PolyhedronTerminalNode : public PolyhedronNode
{
public:
    bool expanded;
    int terminalConstrainIdx;
    PolyhedronTerminalNode(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const int &terminalConstrainIdx);
    bool generateRandomSeed(Eigen::VectorXd& seed, const double &tol, const int &maxTrials);
    void* getNodeData();
private:
    Eigen::VectorXd ai_col;
    double aipivot;
    double bi;
    Eigen::MatrixXd Aj;
    Eigen::VectorXd bj;
    int n;
    int m;
    int col;
    std::vector<int> rows;
    std::vector<int> cols;

    //Eigen::ArithmeticSequence<Eigen::Index,Eigen::Index,Eigen::internal::FixedInt<1>> idxI;
    //Eigen::ArithmeticSequence<Eigen::Index, Eigen::Index, Eigen::internal::FixedInt<1>> idxJ;
};
#endif