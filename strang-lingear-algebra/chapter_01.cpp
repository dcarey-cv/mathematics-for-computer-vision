/*
MATRICES AND GAUSSIAN ELIMINATION
  Learning Objectives:
    - FIRST: implement algorithms and operations manually
    - THEN: compare the results, syntax, and complexity to that of the Eigen library
*/

#include <Eigen/Dense>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <vector>

Eigen::IOFormat eigFmt(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]", "", "", ' ');

void gaussianElimination(Eigen::MatrixXd matrix);

int main()
{
  // Example from Strang Chapter 1.3, page 13
  Eigen::Matrix<double, 3, 4> matrix;
  matrix << 2,  1,  1,  5,
            4, -6,  0, -2,
           -2,  7,  2,  9;
  std::vector<char> variables;

  // My algorithm
  gaussianElimination(matrix);

  // vs. `Eigen::` solution
  Eigen::Matrix3d A = matrix.leftCols(matrix.rows() - 1); // Pulls LHS matrix from my oiginal matrix
  Eigen::Vector3d b = matrix.col(matrix.rows() - 1); // Pulls RHS solution vector from my original matrix
  Eigen::Vector3d eigenSolution = A.lu().solve(b); // Computes the LU decomposition then solves `Ax = b`
  // Alternatively, the above could be written as one line: Eigen::Vector3d eigenSolution = matrix.leftCols(3).lu().solve(matrix.col(3));

  fmt::print("\nComparing...\nEigen library solution:\nu={}, v={}, w={}\n\n", eigenSolution(0), eigenSolution(1), eigenSolution(2));

  return 0;
}

// Precondition: matrix rows must be pivot-ordered
// (no pivot appears beneath a zero in its column)
// since no internal row sorting is performed.
void gaussianElimination(Eigen::MatrixXd matrix)
{
  fmt::print("\nBefore elimination:\n{}\nNote: final column represents RHS\n\n", fmt::streamed(matrix.format(eigFmt)));

  int n = matrix.rows();

  // Forward elimination produces the upper-triangular form of `matrix`
  fmt::print("Eliminating...\n{}\n\t|\n\tV\n", fmt::streamed(matrix.format(eigFmt)));
  for (int col = 0; col < n - 1; ++col)
  {
    for (int row = col + 1; row < n; ++row)
    {
      double factor = matrix(row, col) / matrix(col, col);
      matrix.row(row) -= factor * matrix.row(col);
    }

    fmt::print("{}", fmt::streamed(matrix.format(eigFmt)));
    if(col < (n - 2))
    {
      fmt::print("\n\t|\n\tV\n");
    }
  }

  // Back-substitution solves for w, then v, then u
  // by walking backwards through the reduced matrix
  Eigen::VectorXd solution(n);
  for (int i = n - 1; i >= 0; --i)
  {
    solution(i) = matrix(i, n);
    for (int j = i + 1; j < n; ++j) {
      solution(i) -= matrix(i, j) * solution(j);
    }
    solution(i) /= matrix(i, i);
  }

  fmt::print("\n\nBack-substituting...\nu={}, v={}, w={}\n", solution(0), solution(1), solution(2));
}