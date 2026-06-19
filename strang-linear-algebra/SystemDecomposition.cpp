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
#include <vector>
#include <stdexcept>
#include <numeric>

using fmt::print, fmt::streamed,
      Eigen::MatrixXd, Eigen::VectorXd, Eigen::Index,
      std::vector, std::string, std::abs;

Eigen::IOFormat eigFmt(Eigen::StreamPrecision, 0, ", ", "\n", "[", "]", "", "", ' ');

template<typename Derived>
void threshold(Eigen::MatrixBase<Derived>& mat, double epsilon)
{
  mat = (mat.array().abs() < epsilon).select(0.0, mat);
}

// Ax = b
constexpr double tolerance = 1e-10;
struct SystemOfEquations
{
  MatrixXd CoefficientMatrix; // A
  vector<string> SolutionVector; // x
  VectorXd ConstantVector; // b
  double Epsilon;

  SystemOfEquations(const MatrixXd& A, const vector<string>& vars, const VectorXd& b)
    : CoefficientMatrix(A), SolutionVector(vars), ConstantVector(b),
      Epsilon(tolerance * A.cwiseAbs().maxCoeff())
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
  MatrixXd UpperTriangularMatrix; // U
  VectorXd ConstantVector; // c
  MatrixXd LowerTriangularMatrix; // L
  vector<Index> RowOrder; // Tracks row exchanges
  vector<Index> PivotCols;
  MatrixXd ColumnSpace; // C(A)
  vector<Index> FreeCols;
  unsigned int Rank;
};
RowEchelon toRowEchelon(const SystemOfEquations& sys);
bool pivot(MatrixXd& matrix, vector<Index>& rowOrder, Index col, Index pivotRow, double epsilon, MatrixXd& L);
void eliminate(MatrixXd& matrix, Index col, Index pivotRow, MatrixXd& L);
void buildColumnSpace(const SystemOfEquations& sys, RowEchelon& ref, double epsilon);

// Rx = d
struct ReducedRowEchelon
{ 
  MatrixXd ReducedMatrix; // R
  VectorXd ConstantVector;  // d
  MatrixXd RowSpace; // C(A^T) = C(U^T) = C(R^T)
  MatrixXd Kernel; // N(A)
};
ReducedRowEchelon toReducedRowEchelon(const RowEchelon& ref, double epsilon);
void buildKernel(const RowEchelon& ref, ReducedRowEchelon& rref, double epsilon);


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

/* Page 117: Basic Example
  MatrixXd A(3, 4);
  A <<  1,  3,  3,  2,
        2,  6,  9,  7,
       -1, -3,  3,  4;
  vector<string> vars{"a", "b", "c", "d"};
  VectorXd b(3); b << 0, 0, 0;
*/

///* Page 129: Network Theory
  MatrixXd A(5, 4);
  A << -1,  1,  0,  0,
       -1,  0,  1,  0,
        0, -1,  1,  0,
        0, -1,  0,  1,
        0,  0, -1,  1;
  vector<string> vars{"a", "b", "c", "d"};
  VectorXd b(5); b << 0, 0, 0, 0, 0;
//*/


  SystemOfEquations sys(A, vars, b);
  print("\n\nA =\n{}\n\nx =\n", streamed(sys.CoefficientMatrix.format(eigFmt)));
  for (const string& var : sys.SolutionVector) { print("[{}]\n", var); }
  print("\nb =\n{}\n", streamed(sys.ConstantVector.format(eigFmt)));

  RowEchelon ref = toRowEchelon(sys);
  if (!ref.FreeCols.empty()) { print("\nThe coefficient matrix is rank-deficient.\n"); }
  print("The pivot columns are {}.\n", ref.PivotCols);
  if (std::any_of(ref.RowOrder.begin(), ref.RowOrder.end(), [i = Index(0)](Index j) mutable {return j != i++;}))
  {
    print("Row ordering post-permutation:\n");
    Index middleRow = static_cast<Index>(ref.RowOrder.size() / 2);
    for (Index idx = 0; idx < static_cast<Index>(ref.RowOrder.size()); ++idx)
    {
      print("[row {}] {} [row {}]\n",
        idx,
        (idx == middleRow) ? "->" : "  ",
        ref.RowOrder[idx]);
    }
  }
  buildColumnSpace(sys, ref, sys.Epsilon);
  print("\nU =\n{}\n\nc =\n{}\n\nL =\n{}\n\nC(A) =\n{}\n",
    streamed(ref.UpperTriangularMatrix.format(eigFmt)),
    streamed(ref.ConstantVector.format(eigFmt)),
    streamed(ref.LowerTriangularMatrix.format(eigFmt)),
    streamed(ref.ColumnSpace.format(eigFmt)));

  ReducedRowEchelon rref = toReducedRowEchelon(ref, sys.Epsilon);
  print("\nR =\n{}\n\nd =\n{}\n\nC(A') =\n{}\n\nN(A) =\n{}\n",
    streamed(rref.ReducedMatrix.format(eigFmt)),
    streamed(rref.ConstantVector.format(eigFmt)),
    streamed(rref.RowSpace.format(eigFmt)),
    streamed(rref.Kernel.format(eigFmt)));

  return 0;
}

RowEchelon toRowEchelon(const SystemOfEquations& sys)
{
  RowEchelon ref;

  Index m = sys.CoefficientMatrix.rows();
  Index n = sys.CoefficientMatrix.cols();

  MatrixXd matrix(m, n + 1);
  matrix << sys.CoefficientMatrix, sys.ConstantVector;

  MatrixXd L = MatrixXd::Identity(m, m); // L is always an mxm square matrix
  ref.RowOrder.resize(static_cast<size_t>(m));
  std::iota(ref.RowOrder.begin(), ref.RowOrder.end(), 0);
  ref.PivotCols.reserve(static_cast<size_t>(std::min(m, n)));

  Index pivotRow = 0;
  for (Index col = 0; col < n; ++col)
  {
    if (pivot(matrix, ref.RowOrder, col, pivotRow, sys.Epsilon, L))
    {
      ref.PivotCols.push_back(col);
    } else {
      ref.FreeCols.push_back(col);
      continue;
    }
    
    eliminate(matrix, col, pivotRow++, L);
  }

  ref.Rank = static_cast<unsigned int>(ref.PivotCols.size());
  ref.LowerTriangularMatrix = L;
  ref.UpperTriangularMatrix = matrix.leftCols(n); // U
  ref.ConstantVector = matrix.col(n); // c

  // Thresholding as per the established epsilon
  threshold(ref.UpperTriangularMatrix, sys.Epsilon);
  threshold(ref.ConstantVector, sys.Epsilon);

  return ref; // `ref` stands for "Row Echelon Form"
}

bool pivot(MatrixXd& matrix, vector<Index>& rowOrder, Index col, Index pivotRow, double epsilon, MatrixXd& L)
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

void buildColumnSpace(const SystemOfEquations& sys, RowEchelon& ref, double epsilon)
{
  Index m = sys.CoefficientMatrix.rows();
  Index r = static_cast<Index>(ref.Rank);

  ref.ColumnSpace = MatrixXd::Zero(m, r);
  for (Index idx = 0; idx < r; ++idx) { ref.ColumnSpace.col(idx) = sys.CoefficientMatrix.col(ref.PivotCols[idx]); }

  threshold(ref.ColumnSpace, epsilon);
}

ReducedRowEchelon toReducedRowEchelon(const RowEchelon& ref, double epsilon)
{
  ReducedRowEchelon rref;
  
  Index m = ref.UpperTriangularMatrix.rows();
  Index n = ref.UpperTriangularMatrix.cols();

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

  threshold(rref.ReducedMatrix, epsilon);
  threshold(rref.ConstantVector, epsilon);

  rref.RowSpace = MatrixXd::Zero(n, ref.Rank);
  for(Index idx = 0; idx < ref.Rank; ++idx) { rref.RowSpace.col(idx) = rref.ReducedMatrix.row(idx).transpose(); }
  buildKernel(ref, rref, epsilon);

  return rref; // `rref` stands for "Reduced Row Echelon Form"
}

void buildKernel(const RowEchelon& ref, ReducedRowEchelon& rref, double epsilon)
{
  rref.Kernel = MatrixXd::Zero(rref.ReducedMatrix.cols(), static_cast<Index>(ref.FreeCols.size()));
  for (Index i = 0; i < static_cast<Index>(ref.FreeCols.size()); ++i)
  {
    Index freeCol = ref.FreeCols[i];
    rref.Kernel(freeCol, i) = 1.0;
    for (Index j = 0; j < static_cast<Index>(ref.PivotCols.size()); ++j)
    {
      rref.Kernel(ref.PivotCols[j], i) = -rref.ReducedMatrix(j, freeCol);
    }
  }

  threshold(rref.Kernel, epsilon);
}
