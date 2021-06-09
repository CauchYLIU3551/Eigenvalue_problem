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

std::vector<std::vector<double>> EigenSolver::get_AX(std::vector<std::vector<double>> X )
{
  std::vector<std::vector<double>> AX;
  AX.clear();
  //MX.resize(X.size());
  std::vector<double> temp(X.size(),0);
  for(int j=0;j<X[0].size();j++)
    {
      for (int i=0;i<X.size();i++)
        {
          temp[i]=X[i][j];
        }
      temp=multiply(A,temp);
      AX.push_back(temp);// so in this way, MX[i][j] = MX(j, i) in fact. Because in computing, I store the columns of M*X in every row of MX;
    }
  return AX;
};


// The function is used to compute the projection of v into u 
// under M-inner-productsComputing the project of <u, v> corresponding to M;
std::vector<double> EigenSolver::proj_M(std::vector<double> u, std::vector<double> v)
{
  double delta=0;
  std::vector<double> Mv,Mu;
  Mv=multiply(M,v);
  Mu=multiply(M,u);
  delta=inner(u,Mv)/inner(u,Mu);
  for(int i=0;i<u.size();i++)
    {
      v[i]=delta*u[i];
    }
  return v;
}

// GS orthogonalize the Matrix X corresponding with matrix M;
void EigenSolver::GS_M(std::vector<std::vector<double>> & X)
{
  transpose(X);
  std::vector<std::vector<double>> tempX(X);
  for(int i=1;i<tempX.size();i++)
    {
      for(int j=0;j<i;j++)
        {
          std::vector<double> temp_proj;
          temp_proj=proj_M(X[j], tempX[i]);
          for(int k=0;k<tempX[0].size();k++)
            {
              X[i][k]=X[i][k]-temp_proj[k];
            }
        }
    }

  // std::cout<<"This is the matrix X during the GS_M process\n";
  // show_matrix(X);

  for (int k=0;k<X.size();k++)
    {
      std::vector<double> temp_Xk;
      temp_Xk=multiply(M,X[k]);
      double norm=0;
      norm=sqrt(inner(X[k],temp_Xk));
      for(int i=0;i<X[0].size();i++)
        {
          X[k][i]=X[k][i]/norm;
        }
    }
  transpose(X);
}


void EigenSolver::BIPowerSolve(std::vector<std::vector<double>>& X,
			       std::vector<double>& lambda,
			       int p,
			       int max_iter,
			       double tol)
{
  int k=0;
  double res=10;
  std::vector<double> tempX(A->m(), 0), tempaX(A->m(), 0);
  std::vector<std::vector<double>> H;
  while(k < max_iter)
    {
      transpose(X);
      std::cout<<"The matrix X is :\n";
      show_matrix(X);
      for (int i = 0; i < p; i ++)
	{
	  std::vector<double> tempMb;
	  tempMb = multiply(M, X[i]);
	  CGSolver sol(*A);
	  sol.solve(X[i], tempMb, 1.0e-5, A->m());
	}
      transpose(X);
      std::cout<<"After CGSolver, the result of A^-1 * X is \n";
      show_matrix(X);
      GS_M(X);
      std::cout<<"After M-GS, the result of X is\n";
      show_matrix(X);
      H=get_AX(X);
      std::cout<<"The matrix M*X is \n";
      show_matrix(H);
      multiply(H, X); // Compute the matrix H = Q* A Q;
      std::cout<<"Matrix X' * M *X is \n";
      show_matrix(H);
      QRSolver(H, X);
      std::cout<<"After shur decomposition the matrix H and X are:\n";
      show_matrix(H);
      std::cout<<"Matrix X :\n";
      show_matrix(X);
      // get_residual ;
      for (int i = 0; i < p; i ++)
	{
	  for(int j = 0; j < X.size();j++)
	    {
	      tempX[j] = X[j][i];
	    }
	  tempaX=tempX;
	  tempX=multiply(A, tempX);
	  std::cout<<"The vector AX[i] is \n";
	  show_vector(tempX);
	  tempaX=multiply(M, tempaX);
	  
	  AX(H[i][i], tempaX);
	  std::cout<<"The vector lambda*MX[i] is:\n";
	  show_vector(tempaX);
	  AYPX(-1.0, tempaX, tempX);
	  std::cout<<"After AYPX() the vector AX[i] is :\n";
	  show_vector(tempX);
	  std::cout<<"the lambda*MX[i] is :\n";
	  show_vector(tempaX);
	  std::cout<<"The infi norm is : "<<infi_Norm(tempX)<<"\n";
	  if ( infi_Norm(tempX)<res)
	    {
	      res=infi_Norm(tempX);
	    }
	}
     
      if (res<tol)
	{
	  lambda.clear();
	  for(int j=0; j<H.size();j++)
	    {
	      lambda.push_back(H[j][j]);
	    }
	  std::cout<<"The residual is "<<res<<"\n";
	  return;
	  //break;
	}
      k++;
    }
  lambda.clear();
  for(int j=0; j<H.size();j++)
    {
      lambda.push_back(H[j][j]);
    }
  std::cout<<"The maximum number of iteration exceeded.\n";
  // std::cout<<"The numerical value of eigenvalues are :\n";
}


