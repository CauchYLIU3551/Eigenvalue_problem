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
#include <lac/vector.h>

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
  virtual void get_Ap(std::vector<double>x);
  virtual void get_Ap(dealii::Vector<double>x);
  // compute the residual of the equation Ax=r i.e. res= r-Ax;
  virtual void get_res(std::vector<double> x, std::vector<double> r);
  virtual void get_res(const dealii::Vector<double> x, const dealii::Vector<double> r);
  double tolerence()const {return toler;};

  // solve the equation: Ax=r, to get the x which can satisfy the tolerence;
  void solve(std::vector<double>& x,
	     const std::vector<double> r,
	     double tol=0.0,
	     int max_iter=10);

  void solve(dealii::Vector<double>& x,
	     const dealii::Vector<double>& r,
	     double tol=0.0,
	     int max_iter=10);

  const Matrix* A;
  //std::vector<double> Ap;
  //std::vector<double> res;
  dealii::Vector<double> Ap;
  dealii::Vector<double> res;
 private:
  double toler;
  bool is_initialized;
};


#endif 
