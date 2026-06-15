/*
ELIMINATION AND VECTOR SPACES
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
#include <stdexcept>

using fmt::print,
      Eigen::MatrixXd, Eigen::VectorXd,
      std::vector, std::string;

Eigen::IOFormat eigFmt(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]", "", "", ' ');

// Ax = b
struct SystemOfEquations
{
  MatrixXd CoefficientMatrix; // A
  vector<string> SolutionVector; // x
  VectorXd ConstantVector; // b

  SystemOfEquations(const MatrixXd& A, const vector<string>& vars, const VectorXd& b)
    : CoefficientMatrix(A), SolutionVector(vars), ConstantVector(b)
  {
    if (A.cols() != static_cast<Eigen::Index>(vars.size()))
    // *.cols() returns Eigen::Index (signed ptrdiff_t), whereas *.size() returns size_t (unsigned)
      throw std::invalid_argument("Variable count must match column count of A");
    if (A.rows() != b.size())
      throw std::invalid_argument("b size must match row count of A");
  }
};


// Ux = c
struct RowEchelon
{ 
  MatrixXd UpperTriangularMatrix; // U
  vector<string> SolutionVector; // x
  VectorXd ConstantVector; // c
  vector<int> RowOrder; // Tracks row exchanges
};
RowEchelon toRowEchelon(const SystemOfEquations& sys);

// Rx = d
struct ReducedRowEchelon
{ 
  MatrixXd ReducedMatrix; // R
  vector<string> SolutionVector; // x
  VectorXd ConstantVector;  // d
  vector<int> RowOrder; // Tracks row exchanges
};
ReducedRowEchelon toReducedRowEchelon(const RowEchelon& ref);


// OBJECTIVES:
// - Improve algorithms from chapter_01.cpp in terms of robustness and readability
// - Implement the various stages of decomposition: Ax=b -> Ux=c -> Rx=d
//   - Deduce pivot variables and free variables
//     - Determine the rank of the coefficient matrix, A, via its reductions, U and R
//     - Extract Nullspace solution vectors and particular solutions
int main()
{
  MatrixXd A(3, 4); // SYSTEM OF EQUATIONS
  A << 1, 2, 3, 5,  // 1x1 + 2x2 + 3x3 + 5x4 = b1;
       2, 4, 8, 12, // 2x1 + 4x2 + 8x3 + 12x4 = b2;
       3, 6, 7, 13; // 3x1 + 6x2 + 7x3 + 13x4 = b3;

  vector<string> vars{"x1", "x2", "x3", "x4"};

  VectorXd b(3);
  b << 0, 6, -6;  // b1, b2, b3

  SystemOfEquations sys(A, vars, b);
  RowEchelon ref = toRowEchelon(sys);
  ReducedRowEchelon rref = toReducedRowEchelon(ref);

  return 0;
}

RowEchelon toRowEchelon(const SystemOfEquations& sys)
{
  RowEchelon ref;

  return ref;
}

ReducedRowEchelon toReducedRowEchelon(const RowEchelon& ref)
{
  ReducedRowEchelon rref;

  return  rref;
}