#ifndef _EigenSolver_h_
#define _EigenSolver_h_

#include <lac/sparse_matrix.h>
#include <lac/sparsity_pattern.h>
#include <math.h>
#include <vector>
//#include<EigenSolver/Miscellaneous.h>
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

  // BIPowerSolve(): Block Inverse Power Methods to compute the p smallest eigenpairs. x represents the eigenvectors,
  // p represents the number of eigenpairs,
  // lambda represents the p smallest eigenvalues.
  void BIPowerSolve(std::vector<std::vector<double>>& X,
		    std::vector<double>&lambda,
		    int p=2,
		    int max_iter = 100,
		    double tol = 1.e-5);

  // LOBPCGSolve(): LOBPCG Methods to compute p smallest eigenpairs-
  // of matrix A corresponding to M.
  void LOBPCGSolve(std::vector<std::vector<double>>& X,
		   std::vector<double>&lambda,
		   int p=2,
		   int max_iter = 100,
		   double tol = 1.e-5);

  // Get (AX)T;
  std::vector<std::vector<double>> get_AX(std::vector<std::vector<double>> X);
  // Get (MX)T;
  std::vector<std::vector<double>> get_MX(std::vector<std::vector<double>> X);

  // Get M^-1 * A * X;
  std::vector<std::vector<double>> get_MAX(std::vector<std::vector<double>> X);
  // The function is used to compute the projection of v into u 
  // under M-inner-productsComputing the project of <u, v> corresponding to M;
  std::vector<double> proj_M(std::vector<double> u, std::vector<double> v);

  // The function is used to compute the projection of v into u 
  // under A-inner-productsComputing the project of <u, v> corresponding to A;
  std::vector<double> proj_A(std::vector<double> u, std::vector<double> v);
  
  // GS orthogonalize the Matrix X corresponding with matrix M;
  void GS_M(std::vector<std::vector<double>> & X);
  // GS orthogonalize the Matrix X corresponding with matrix A;
  void GS_A(std::vector<std::vector<double>> & X);
 private:
  const Matrix* A;
  const Matrix* M;
  
};

#endif
