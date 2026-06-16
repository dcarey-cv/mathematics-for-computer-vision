/*
ELIMINATION AND VECTOR SPACES
  Learning Objectives:
    - FIRST: implement algorithms and operations manually
    - THEN: compare the results, syntax, and complexity to that of the Eigen library
*/

#include <Eigen/Dense>

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <string>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <numeric>

using fmt::print, fmt::streamed,
      Eigen::MatrixXd, Eigen::VectorXd, Eigen::Index,
      std::vector, std::string, std::abs;

Eigen::IOFormat eigFmt(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]", "", "", ' ');

// Ax = b
constexpr double tolerance = 1e-10;
struct SystemOfEquations
{
  MatrixXd CoefficientMatrix; // A
  vector<string> SolutionVector; // x
  VectorXd ConstantVector; // b
  unsigned int Rank; // r

  SystemOfEquations(const MatrixXd& A, const vector<string>& vars, const VectorXd& b)
    : CoefficientMatrix(A), SolutionVector(vars), ConstantVector(b)
  {
    if (A.cols() != static_cast<Index>(vars.size()))
      // *.cols() returns Eigen::Index (signed ptrdiff_t),
      //whereas *.size() returns size_t (unsigned)
      throw std::invalid_argument("Variable count must match column count of A");
    if (A.rows() != b.size())
      throw std::invalid_argument("b size must match row count of A");
  }
};

// Ux = c
struct RowEchelon
{ 
  MatrixXd LowerTriangularMatrix; // L
  MatrixXd UpperTriangularMatrix; // U
  VectorXd ConstantVector; // c
  bool IsPermuted = false;
  vector<Index> RowOrder; // Tracks row exchanges
  vector<Index> PivotCols;
  unsigned int Rank;
};
RowEchelon toRowEchelon(const SystemOfEquations& sys);
bool pivot(MatrixXd& matrix, vector<Index>& rowOrder, Index col, Index pivotRow, double epsilon, MatrixXd& L, bool& isPermuted);
void eliminate(MatrixXd& matrix, Index col, Index pivotRow, MatrixXd& L);

// Rx = d
struct ReducedRowEchelon
{ 
  MatrixXd ReducedMatrix; // R
  VectorXd ConstantVector;  // d
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
/* Page 94: Another Worked Example
  MatrixXd A(3, 4); // SYSTEM OF EQUATIONS
  A << 1, 2, 3, 5,  // 1x1 + 2x2 + 3x3 + 5x4 = b1;
       2, 4, 8, 12, // 2x1 + 4x2 + 8x3 + 12x4 = b2;
       3, 6, 7, 13; // 3x1 + 6x2 + 7x3 + 13x4 = b3;
  vector<string> vars{"x1", "x2", "x3", "x4"};
  VectorXd b(3); b << 0, 6, -6;  // b1, b2, b3
*/

///* Problem Set 2.2: 2A
  MatrixXd A(2, 4);
  A << 0, 1, 0, 3,
       0, 2, 0, 6;
  vector<string> vars{"a", "b", "c", "d"};
  VectorXd b(2); b << 0, 0;  // b1, b2
//*/


  SystemOfEquations sys(A, vars, b);
  print("\n\nA =\n{}\n\nx =\n", streamed(sys.CoefficientMatrix.format(eigFmt)));
  for (const string& var : sys.SolutionVector) { print("[{}]\n", var); }
  print("\nb =\n{}\n", streamed(sys.ConstantVector.format(eigFmt)));

  RowEchelon ref = toRowEchelon(sys);
  sys.Rank = ref.Rank;
  if (sys.Rank < sys.SolutionVector.size()) { print("\nThe coefficient matrix is rank-deficient.\n"); }
  print("The pivot columns are {}.\n", ref.PivotCols);
  if (ref.IsPermuted)
  {
    print("Row ordering post-permutation:\n");
    for (Index idx = 0; idx < static_cast<Index>(ref.RowOrder.size()); ++idx)
    {
      print("[row {}] {} [row {}]\n",
        idx,
        (idx == static_cast<Index>(ref.RowOrder.size() / 2)) ? "->" : "  ",
        ref.RowOrder[idx]);
    }
  }
  print("\nU =\n{}\n\nc =\n{}\n\nL =\n{}\n",
    streamed(ref.UpperTriangularMatrix.format(eigFmt)),
    streamed(ref.ConstantVector.format(eigFmt)),
    streamed(ref.LowerTriangularMatrix.format(eigFmt)));

  ReducedRowEchelon rref = toReducedRowEchelon(ref);
  print("\nR =\n{}\n\nd =\n{}\n\n",
    streamed(rref.ReducedMatrix.format(eigFmt)),
    streamed(rref.ConstantVector.format(eigFmt)));

  return 0;
}



RowEchelon toRowEchelon(const SystemOfEquations& sys)
{
  RowEchelon ref;

  Index m = sys.CoefficientMatrix.rows();
  Index n = sys.CoefficientMatrix.cols();
  double epsilon = tolerance * sys.CoefficientMatrix.cwiseAbs().maxCoeff();
    // Relative epsilon based on largest coefficient

  MatrixXd matrix(m, n + 1);
  matrix << sys.CoefficientMatrix, sys.ConstantVector;

  MatrixXd L = MatrixXd::Identity(m, m); // L is always an mxm square matrix
  ref.RowOrder.resize(static_cast<size_t>(m));
  std::iota(ref.RowOrder.begin(), ref.RowOrder.end(), 0);
  ref.PivotCols.reserve(static_cast<size_t>(std::min(m, n)));

  Index pivotRow = 0;
  for (Index col = 0; col < n; ++col)
  {
    if (!pivot(matrix, ref.RowOrder, col, pivotRow, epsilon, L, ref.IsPermuted)) { continue; }
    ref.PivotCols.push_back(col);
    eliminate(matrix, col, pivotRow++, L);
  }

  ref.Rank = static_cast<unsigned int>(ref.PivotCols.size());
  ref.LowerTriangularMatrix = L;
  ref.UpperTriangularMatrix = matrix.leftCols(n); // U
  ref.ConstantVector = matrix.col(n); // c

  // Thresholding as per the established epsilon
  ref.UpperTriangularMatrix = (ref.UpperTriangularMatrix.array().abs() < epsilon).select(0.0, ref.UpperTriangularMatrix);
  ref.ConstantVector = (ref.ConstantVector.array().abs() < epsilon).select(0.0, ref.ConstantVector);

  return ref; // `ref` stands for "Row Echelon Form"
}

bool pivot(MatrixXd& matrix, vector<Index>& rowOrder, Index col, Index pivotRow, double epsilon, MatrixXd& L, bool& isPermuted)
{
  Index maxRow = pivotRow;
  for (Index idx = maxRow + 1; idx < matrix.rows(); ++idx)
  {
    if (abs(matrix(idx, col)) <= abs(matrix(maxRow, col))) { continue; }
    else { maxRow = idx; }
  }

  if (abs(matrix(maxRow, col)) < epsilon) { return false; }

  if (maxRow != pivotRow)
  {
    isPermuted = true;
    matrix.row(pivotRow).swap(matrix.row(maxRow));
    L.block(pivotRow, 0, 1, pivotRow).swap(L.block(maxRow, 0, 1, pivotRow));
    std::swap(rowOrder[pivotRow], rowOrder[maxRow]);
  }

  return true;
}

void eliminate(MatrixXd& matrix, Index col, Index pivotRow, MatrixXd& L)
{
  for (Index idx = pivotRow + 1; idx < matrix.rows(); ++idx)
  {
    double factor = matrix(idx, col) / matrix(pivotRow, col);
    matrix.row(idx) -= factor * matrix.row(pivotRow);
    L(idx, pivotRow) =  factor;
  }
}



ReducedRowEchelon toReducedRowEchelon(const RowEchelon& ref)
{
  ReducedRowEchelon rref;
  
  Index m = ref.UpperTriangularMatrix.rows();
  Index n = ref.UpperTriangularMatrix.cols();
  double epsilon = tolerance * ref.UpperTriangularMatrix.cwiseAbs().maxCoeff();

  MatrixXd matrix(m, n + 1);
  matrix << ref.UpperTriangularMatrix, ref.ConstantVector;

  for (Index pivotRow = ref.Rank - 1; pivotRow >= 0; --pivotRow)
  {
    Index pivotCol = ref.PivotCols[pivotRow];
    matrix.row(pivotRow) /= matrix(pivotRow, pivotCol);
    for (Index j = 0; j < pivotRow; ++j)
    {
      matrix.row(j) -= matrix(j, pivotCol) * matrix.row(pivotRow);
    }
  }

  rref.ReducedMatrix = matrix.leftCols(n);
  rref.ConstantVector = matrix.col(n);

  rref.ReducedMatrix = (rref.ReducedMatrix.array().abs() < epsilon).select(0.0, rref.ReducedMatrix);
  rref.ConstantVector = (rref.ConstantVector.array().abs() < epsilon).select(0.0, rref.ConstantVector);

  return rref; // `rref` stands for "Reduced Row Echelon Form"
}