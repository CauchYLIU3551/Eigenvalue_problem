#ifndef _CGSolver_h_
#define _CGSolver_h_

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <assert.h> 

#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>

class CGSolver
{
 public:
  typedef dealii::SparseMatrix<double> Matrix;
  CGSolver();
  CGSolver(const Matrix&,
  	   double tol = 1.0e-12);
  virtual ~CGSolver();

  void reinit(const Matrix&);
  void clear();

  // compute the value of matrix A times p;
  void get_Ap(std::vector<double>x);

  // compute the residual of the equation Ax=r i.e. res= r-Ax;
  void get_res(const std::vector<double> x, const std::vector<double> r);

  double tolerence()const {return toler;};

  // solve the equation: Ax=r, to get the x which can satisfy the tolerence;
  void solve(std::vector<double>& x,
	     const std::vector<double> r,
	     double tol=0.0,
	     int max_iter=10);

  const Matrix* A;
  std::vector<double> Ap;
  std::vector<double> res;
 private:
  double toler;
  bool is_initialized;
};


#endif 
