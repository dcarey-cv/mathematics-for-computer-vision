/*
MATRICES AND GAUSSIAN ELIMINATION
  Learning Objectives:
    - FIRST: implement algorithms and operations manually
    - THEN: compare the results, syntax, and complexity to that of the Eigen library
*/

#include <Eigen/Dense>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <string>
#include <unordered_map>
#include <vector>

Eigen::IOFormat eigFmt(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]", "", "", ' ');

struct EliminationStep { int targetRow, pivotCol; double factor; };
struct PermutationStep { int row1, row2; };
struct ElimResult {
  Eigen::VectorXd solution;
  std::vector<std::string> variables;
  std::vector<EliminationStep> eliminations;
  std::vector<PermutationStep> permutations;
  Eigen::MatrixXd U;
  bool valid;
};

ElimResult gaussianElimination(Eigen::MatrixXd matrix, std::vector<std::string> variables);
EliminationStep eliminate(const Eigen::MatrixXd& matrix, const int& pivotCol, const int& targetRow);
bool partialPivot(Eigen::MatrixXd& matrix, const int& col, std::vector<PermutationStep>& permutations, std::vector<std::string>& variables);
Eigen::MatrixXd reconstructL(const ElimResult& result);
Eigen::MatrixXd reconstructD(const ElimResult& result);
Eigen::MatrixXd reconstructU(const ElimResult& result);
bool gaussElimComparison(const Eigen::MatrixXd& matrix, const ElimResult& result, const std::vector<std::string>& originalVars);

int main()
{
  // Example from Strang Chapter 1.3, page 13
  Eigen::Matrix<double, 3, 4> matrix;
  matrix << 2,  1,  1,  5,
            4, -6,  0, -2,
           -2,  7,  2,  9;
  std::vector<std::string> variables;
  variables.push_back("u");
  variables.push_back("v");
  variables.push_back("w");

  int bCol = matrix.cols() - 1;
  int midRow = matrix.rows() / 2;
  for (int idx = 0; idx < matrix.rows(); ++idx)
  {
    fmt::print("[");
    for (int j = 0; j < bCol; ++j)
    {
      fmt::print("{:4g}{}", matrix(idx, j), (j < bCol - 1) ? ", " : "");
    }
    fmt::print("] [{:>2}] {} [{:4g}]\n", variables[idx], (idx == midRow) ? "=" : " ", matrix(idx, bCol));
  }

  // My algorithm
  ElimResult result = gaussianElimination(matrix, variables);
  if (!result.valid) { return 1; }

  std::unordered_map<std::string, double> solutionMap;
  for (int idx = 0; idx < static_cast<int>(result.variables.size()); ++idx)
  {
    solutionMap[result.variables[idx]] = result.solution(idx);
  }
  for (const std::string& var : variables)
  {
    fmt::print("{}={}\t", var, solutionMap[var]);
  }
  gaussElimComparison(matrix, result, variables); // Verify against Eigen library solution

  // LDU decomposition
  fmt::print("\nL:\n{}\n\nD:\n{}\n\nU:\n{}\n",
             fmt::streamed(reconstructL(result).format(eigFmt)),
             fmt::streamed(reconstructD(result).format(eigFmt)),
             fmt::streamed(reconstructU(result).format(eigFmt)));

  return 0;
}


ElimResult gaussianElimination(Eigen::MatrixXd matrix, std::vector<std::string> variables)
{
  int n = matrix.rows();

  if (static_cast<int>(variables.size()) != n)
  {
    fmt::print("Error: {} variable names provided but matrix has {} rows\n", variables.size(), n);
    return { {}, {}, {}, {}, {}, false };
  }

  std::vector<EliminationStep> eliminations;
  std::vector<PermutationStep> permutations;

  for (int col = 0; col < n - 1; ++col)
  {
    partialPivot(matrix, col, permutations, variables);

    for (int row = col + 1; row < n; ++row)
    {
      EliminationStep step = eliminate(matrix, col, row);
      matrix.row(row) -= step.factor * matrix.row(col);
      eliminations.push_back(std::move(step));
    }
  }

  Eigen::VectorXd solution(n);
  for (int idx = n - 1; idx >= 0; --idx)
  {
    solution(idx) = matrix(idx, n);
    for (int j = idx + 1; j < n; ++j)
    {
      solution(idx) -= matrix(idx, j) * solution(j);
    }
    solution(idx) /= matrix(idx, idx);
  }

  return { std::move(solution), std::move(variables), std::move(eliminations), std::move(permutations), matrix.leftCols(n), true };
}

bool partialPivot(Eigen::MatrixXd& matrix, const int& col, std::vector<PermutationStep>& permutations, std::vector<std::string>& variables)
{
  int n = matrix.rows();
  int maxRow = col;
  const int& pivotRow = col;

  for (int idx = pivotRow + 1; idx < n; ++idx)
  {
    if (std::abs(matrix(idx, pivotRow)) > std::abs(matrix(maxRow, pivotRow)))
    {
      maxRow = idx;
    }
  }

  if (maxRow == pivotRow)
  {
    return false;
  }

  matrix.row(pivotRow).swap(matrix.row(maxRow));
  std::swap(variables[pivotRow], variables[maxRow]);
  permutations.push_back({ pivotRow, maxRow });
  return true;
}

EliminationStep eliminate(const Eigen::MatrixXd& matrix, const int& pivotCol, const int& targetRow)
{
  return { targetRow, pivotCol, matrix(targetRow, pivotCol) / matrix(pivotCol, pivotCol) };
}

Eigen::MatrixXd reconstructL(const ElimResult& result)
{
  int n = result.U.rows();
  Eigen::MatrixXd L = Eigen::MatrixXd::Identity(n, n);
  for (const EliminationStep& step : result.eliminations)
  {
    L(step.targetRow, step.pivotCol) = step.factor;
  }
  return L;
}

Eigen::MatrixXd reconstructD(const ElimResult& result)
{
  int n = result.U.rows();
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(n, n);
  for (int idx = 0; idx < n; ++idx)
  {
    D(idx, idx) = result.U(idx, idx);
  }
  return D;
}

Eigen::MatrixXd reconstructU(const ElimResult& result)
{
  Eigen::MatrixXd Uunit = result.U;
  int n = Uunit.rows();
  for (int idx = 0; idx < n; ++idx)
  {
    Uunit.row(idx) /= result.U(idx, idx);
  }
  return Uunit;
}


bool gaussElimComparison(const Eigen::MatrixXd& matrix, const ElimResult& result, const std::vector<std::string>& originalVars)
{
  int n = matrix.cols() - 1;

  Eigen::VectorXd eigenSolution = matrix.leftCols(n).lu().solve(matrix.col(n));

  std::unordered_map<std::string, double> solutionMap;
  for (int idx = 0; idx < n; ++idx)
  {
    solutionMap[result.variables[idx]] = result.solution(idx);
  }

  const double tol = 1e-9;
  for (int idx = 0; idx < n; ++idx)
  {
    if (std::abs(solutionMap[originalVars[idx]] - eigenSolution(idx)) > tol)
    {
      fmt::print("VERIFICATION FAILED: {}={} (expected {})\n",
                 originalVars[idx], solutionMap[originalVars[idx]], eigenSolution(idx));
      return false;
    }
  }

  fmt::print("\nVerified.\n");
  return true;
}