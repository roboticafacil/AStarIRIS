#include<Eigen/Dense>
#include "EigenUtils.h"
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;

std::shared_ptr<ndarray<double, 2>> Eigen2NdArray(const Eigen::MatrixXd& A)
{
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Arow(A);
    std::shared_ptr<ndarray<double, 2>> A_ptr(new ndarray<double, 2>(Arow.data(), shape(A.rows(), A.cols())));
    return A_ptr;
}

std::shared_ptr<ndarray<double, 1>> Eigen2NdArray(const Eigen::VectorXd& b)
{
    std::shared_ptr<ndarray<double, 1>> b_ptr(new ndarray<double, 1>(b.data(), shape(b.rows())));
    return b_ptr;
}

double nchoosek(const int n, const int k)
{
    if (k > n / 2)
    {
        return nchoosek(n, n - k);
    }
    else if (k == 1)
    {
        return n;
    }
    else
    {
        double c = 1;
        for (int i = 1; i <= k; i++)
        {
            c *= (((double)n - k + i) / ((double)i));
        }
        return std::round(c);
    }
}

void nchoosek(const Eigen::VectorXi& V, const int k, Eigen::MatrixXi& U)
{
    if (V.size() == 0)
    {
        U.resize(0, k);
        return;
    }
    assert((V.cols() == 1 || V.rows() == 1) && "V must be a vector");
    U.resize(nchoosek(V.size(), k), k);
    int running_i = 0;
    int running_j = 0;
    Eigen::VectorXi running(k);
    int N = V.size();
    const std::function<void(int, int)> doCombs =
        [&running, &N, &doCombs, &running_i, &running_j, &U, &V](int offset, int k)
    {
        if (k == 0)
        {
            U.row(running_i) = running;
            running_i++;
            return;
        }
        for (int i = offset; i <= N - k; ++i)
        {
            running(running_j) = V(i);
            running_j++;
            doCombs(i + 1, k - 1);
            running_j--;
        }
    };
    doCombs(0, k);
    U = Eigen::MatrixXi(U.array() - 1);
}

// Function to compare two rows with tolerance due to floating-point precision issues
bool areRowsEqual(const Eigen::RowVectorXd& row1, const Eigen::RowVectorXd& row2, const double& tolerance) {
    return (row1 - row2).norm() < tolerance; // Compare the norm of the difference between rows
}

// Function to remove duplicate rows from an Eigen::MatrixXd
Eigen::MatrixXd removeDuplicateRows(const Eigen::MatrixXd& matrix, const double& tolerance) {
    std::vector<Eigen::RowVectorXd> uniqueRows;

    for (int i = 0; i < matrix.rows(); ++i) {
        Eigen::RowVectorXd currentRow = matrix.row(i);
        bool isDuplicate = false;

        // Check if the current row already exists in uniqueRows
        for (const auto& uniqueRow : uniqueRows) {
            if (areRowsEqual(currentRow, uniqueRow, tolerance)) {
                isDuplicate = true;
                break;
            }
        }

        // If not a duplicate, add it to uniqueRows
        if (!isDuplicate) {
            uniqueRows.push_back(currentRow);
        }
    }

    // Build the new matrix from uniqueRows
    Eigen::MatrixXd result(uniqueRows.size(), matrix.cols());
    for (size_t i = 0; i < uniqueRows.size(); ++i) {
        result.row(i) = uniqueRows[i];
    }
    return result;
}

void unique(Eigen::MatrixXd& A)
{
    std::vector<Eigen::VectorXd> vec;
    for (int64_t i = 0; i < A.rows(); ++i)
        vec.push_back(A.row(i));

    std::sort(vec.begin(), vec.end(), [](Eigen::VectorXd const& t1, Eigen::VectorXd const& t2) { return t1(0) < t2(0); });

    auto it = std::unique(vec.begin(), vec.end());
    vec.resize(std::distance(vec.begin(), it));

    A.resize(vec.size(), A.cols());
    for (int64_t i = 0; i < vec.size(); ++i)
        A.row(i) = vec[i];
}