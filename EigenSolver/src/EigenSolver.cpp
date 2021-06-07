#include<EigenSolver/EigenSolver.h>
#include<EigenSolver/Miscellaneous.h>

#include "CG/CGSolver.h"

#include<stdio.h>
#include<assert.h>
#include<iostream>
EigenSolver::EigenSolver(){};
EigenSolver::EigenSolver(const Matrix& a, const Matrix& m){
  A=&a;
  M=&m;
  assert(A->m()==A->n());
  assert(M->m()==M->n());
  assert(A->m()==M->m());
  // initializing the vector of theta;
  // theta.resize(A->m());
}

void showvector(std::vector<double> x)
{
  for(int i=0;i<x.size();i++)
    {
      std::cout<<x[i]<<" ";
    }
  std::cout<<"\n";
}

void EigenSolver::PowerSolve(std::vector<double>& x,
			     double lambda,
			     int max_iter,
			     double tol)
{
  int iter = 0;
  int p;
  double ERR;
  std::vector<double> y(A->m()),tempAx;
  max_norm(x, p);
  lambda = x[p];
  //std::cout<<"Lambda is "<<lambda<<"\n";
  AX(1.0/lambda, x);
  // showvector(x);
  while (iter < max_iter)
    {
      // std::cout<<"\n\nThis is the "<<iter<<"-th iteration!!!\n";
      tempAx = multiply(A, x);
      //std::cout<<"This is the vector x\n";
      // showvector(x);
      
      CGSolver sol(*M);
      sol.solve(y, tempAx, 1.0e-3,20);

      // std::cout<<"This is vector y\n";
      //showvector(y);
      
      max_norm(y, p);
      lambda = y[p];
      //std::cout<<"Lambda is "<<lambda<<"\n";
      if (fabs(lambda) < tol)
	{
	  std::cout<<"Lambda is "<<lambda<<"\n And the vector y is:\n";
	  std::cout<<"The matrices have the eigenvalue 0.\n";
	  break;
	}
      AX(1.0/lambda, y);
      //std::cout<<"This is vector y/y_p : \n";
      // showvector(y);      
      AYPX(-1.0, y, x);
      // std::cout<<"This is new x\n";
      //showvector(x);
      ERR = infi_Norm(x);
      //std::cout<<"The error is "<<ERR<<"\n";
      x=y;     
      if (ERR < tol)
	{
	  std::cout<<"After "<<iter+1<<" iterations, the eigenvalue is "<<lambda<<"\n";
	  break;
	}
      iter++;
    }
  if (iter == max_iter)
    {
      std::cout<<"The maximum number of iterations exceeded";
    }
}
