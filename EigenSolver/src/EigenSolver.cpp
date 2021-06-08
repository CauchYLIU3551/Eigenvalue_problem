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

//////////////////////////////////////////////
/*  This function is Power method solver to get the maximum eigenvalue and related *  eigenvector of the matrix A corresponding to matrix M. In classic case, 
 *  matrix M can be set as an identity matrix.
 */
void EigenSolver::PowerSolve(std::vector<double>& x,
			     double& lambda,
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
     
      tempAx=multiply(A,x);
      //std::cout<<"This is the vector x\n";
      // showvector(x);
      
      CGSolver sol(*M);
      sol.solve(y, tempAx, 1.0e-5,A->m());

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
	  std::cout<<"After "<<iter+1<<" iterations, the maximum eigenvalue is "<<lambda<<"\n";
	  break;
	}
      iter++;
    }
  if (iter == max_iter)
    {
      std::cout<<"The maximum number of iterations exceeded and the eigenvalue is"<<lambda<<"\n";
    }
}

//////////////////////////////////////////////////////////////////////
/*
 *  This is the Inverse Power Method solver to get the minimum eigenvalue and 
 *  related eigenvector of the Matrix A corresponding to M. In classic case, 
 *  M can be set as identity matrix.
 */
void EigenSolver::IPowerSolve(std::vector<double>& x,
			      double& lambda,
			      int max_iter,
			      double tol)
{
  int iter = 0;
  int p;
  double ERR;
  std::vector<double> y(A->m()),tempMx,tempyyy;
  max_norm(x, p);
  lambda = x[p];
  // std::cout<<"Lambda is "<<lambda<<"\n";
  AX(1.0/lambda, x);
  //showvector(x);
  while (iter < max_iter)
    {
      //  std::cout<<"\n\nThis is the "<<iter<<"-th iteration!!!\n";
      tempMx = multiply(M, x);
      
      // std::cout<<"This is the vector x\n";
      // showvector(x);
      
      CGSolver sol(*A);
      sol.solve(y, tempMx, 1.0e-5,A->m());

      // std::cout<<"This is vector y\n";
      // showvector(y);
      
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
      // std::cout<<"This is vector y/y_p : \n";
      // showvector(y);
      tempyyy=x;
      AYPX(-1.0, y, x);
      // std::cout<<"This is new x\n";
      // showvector(x);
      ERR = infi_Norm(x);
      // std::cout<<"The error is "<<ERR<<"\n";
      x=y;     
      if (ERR < tol)
	{
	  lambda=1.0/lambda;
	  std::cout<<"After "<<iter+1<<" iterations, the minimum eigenvalue is "<<lambda<<"\n";	  
	  break;
	}
      iter++;
    }
  if (iter == max_iter)
    {
      std::cout<<"The eigen is "<<1.0/lambda<<"\n";
      std::cout<<"The ERR is "<<ERR<<"\n";
      std::cout<<"The vector x is \n";
      showvector(tempyyy);
      std::cout<<"The vector y is \n";
      showvector(y);
      std::cout<<"The maximum number of iterations exceeded\n";
    }
}
