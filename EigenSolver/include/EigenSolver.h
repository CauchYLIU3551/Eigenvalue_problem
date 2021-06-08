#ifndef _EigenSolver_h_
#define _EigenSolver_h_

#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <math.h>
#include <vector>
class EigenSolver
{
 public:
  typedef dealii::SparseMatrix<double> Matrix;
  EigenSolver();
  EigenSolver(const Matrix& a, const Matrix& m);
  void PowerSolve(std::vector<double>& x, 
		  double& lambda, 
		  int max_iter = 100,
		  double tol = 1.e-3);
  void IPowerSolve(std::vector<double>& x,
		   double& lambda,
		   int max_iter = 100,
		   double tol = 1.e-3);
 private:
  const Matrix* A;
  const Matrix* M;
  
};

#endif
