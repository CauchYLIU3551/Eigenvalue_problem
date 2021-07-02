#include<EigenSolver/EigenSolver.h>
#include<EigenSolver/Miscellaneous.h>

#include "CG/CGSolver.h"

#include<stdio.h>
#include<assert.h>
#include<iostream>
#include<algorithm>
#include<numeric>

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

std::vector<std::vector<double>> EigenSolver::get_MX(std::vector<std::vector<double>> X )
{
  std::vector<std::vector<double>> MX;
  MX.clear();
  //MX.resize(X.size());
  std::vector<double> temp(X.size(),0);
  for(int j=0;j<X[0].size();j++)
    {
      for (int i=0;i<X.size();i++)
        {
          temp[i]=X[i][j];
        }
      temp=multiply(M,temp);
      MX.push_back(temp);// so in this way, MX[i][j] = MX(j, i) in fact. Because in computing, I store the columns of M*X in every row of MX;
    }
  return MX;
};

// This function compute M^-1 * A * X;
std::vector<std::vector<double>> EigenSolver::get_MAX(std::vector<std::vector<double>> X )
{
  std::vector<std::vector<double>> MAX;
  MAX.clear();
  //MX.resize(X.size());
  std::vector<double> temp(X.size(),0), tempAx;
  CGSolver sol(*M);
  for(int j=0;j<X[0].size();j++)
    {
      for (int i=0;i<X.size();i++)
        {
          temp[i]=X[i][j];
        }
      tempAx=multiply(A,temp);
      sol.solve(temp, tempAx, 1.e-5, A->m()); // This step solve that M temp = tempAx, to get M^-1 * A * x;
      MAX.push_back(temp);// so in this way, MX[i][j] = MX(j, i) in fact. Because in computing, I store the columns of M*X in every row of MX;
    }
  transpose(MAX);
  return MAX;
}

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

// The function is used to compute the projection of v into u 
// under A-inner-productsComputing the project of <u, v> corresponding to A;
std::vector<double> EigenSolver::proj_A(std::vector<double> u, std::vector<double> v)
{
  double delta=0;
  std::vector<double> Av,Au;
  Av=multiply(A,v);
  Au=multiply(A,u);
  delta=inner(u,Av)/inner(u,Au);
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

// GS orthogonalize the Matrix X corresponding with matrix A;
void EigenSolver::GS_A(std::vector<std::vector<double>> & X)
{
  transpose(X);
  std::vector<std::vector<double>> tempX(X);
  for(int i=1;i<tempX.size();i++)
    {
      for(int j=0;j<i;j++)
        {
          std::vector<double> temp_proj;
          temp_proj=proj_A(X[j], tempX[i]);
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
      temp_Xk=multiply(A,X[k]);
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
  double res=10, temp_eig;
  std::vector<double> tempMb, tempR(p);
  std::vector<double> tempX(A->m(), 0), tempaX(A->m(), 0);
  std::vector<std::vector<double>> H, Qk, MAX, tempA, R(p, tempR), Z, tempZ(p,tempR),tempQ(p,tempR);
  GS_M(X);

  while(k < max_iter)
    {
      transpose(X);
      //std::cout<<"The matrix X is :\n";
      //show_matrix(X);
      for (int i = 0; i < p; i ++)
	{
	  tempMb = multiply(M, X[i]);
	  CGSolver sol(*A);
	  sol.solve(X[i], tempMb, 1.0e-5, A->m());
	}
      transpose(X);
      //Z = X;

      GS_M(X);


      /*Following part are the algorithm of subspace iteration method to get several eigenpairs;
       *
       */
            /************************************************
      for (int i=0; i<p; i++)
	{
	  for(int j=0; j<p; j++)
	    {
	      tempZ[i][j] = Z[i][j];
	      tempQ[i][j] = X[i][j];
	    }
	}

      transpose(tempZ);
      for(int i=0;i<p;i++)
	{
	  Gauss(tempQ, R[i], tempZ[i]);
	}
      transpose(tempZ);


      res = 0;
      transpose(X);
      for(int i = 0;i < p ; i ++)
	{
	  tempMb =multiply(M,X[i]);
	  CGSolver sol2(*A);
	  sol2.solve(tempaX, tempMb, 1.0e-5,A->m());
	  temp_eig = R[i][i];
	  AYPX(-1.0*temp_eig, X[i], tempaX);
	  if(res<infi_Norm(tempaX))
	    {
	      res = infi_Norm(tempaX);
	    }
	}

      //  std::cout<<"The residual is ::\n";
      // std::cout<<res<<std::endl;
      transpose(X);

            */////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      // std::cout<<"After M-GS, the result of X is\n";
      // show_matrix(X);


      ////// Begin to compute X'* M^-1 * A *X;
      H=get_AX(X);
      multiply(H,X);
      // std::cout<<"The matrix A*X is \n";
      // show_matrix(H);


      ////// Compute M^-1 * A *X
      
      /////////////// The R-R step in this algorithm might be useless, Ritz vectors converges to eigenvectors but not the p smallest; 
      //////// Then we can get the matrix M^-1 * A * X;
      
      // transpose(X);
      //  H = X;
      // transpose(X);
      // MAX = get_MAX(X);
      
      //multiply(H, MAX); // Compute the matrix H = Q* A Q;
      // std::cout<<"Matrix X' * M *X is \n";
      // show_matrix(H);
      //  std::cout<<"The size of matrix X'*M*X is "<<H.size()<<" x "<<H[0].size()<<"\n";
      identitymatrix(Qk, p);
      // std::cout<<"The matrix Qk is \n";
      // show_matrix(Qk);
      QRSolver(H, Qk);
      // std::cout<<"After shur decomposition the matrix H and X are:\n";
      // show_matrix(H);
      multiply(X,Qk);
      // std::cout<<"Matrix X :\n";
      // show_matrix(X);
      // get_residual ;
     

      
      res=0;
      for (int i = 0; i < p; i ++)
	{
	  for(int j = 0; j < X.size();j++)
	    {
	      tempX[j] = X[j][i];
	    }
	  tempaX=tempX;
	  tempX=multiply(A, tempX);
	  //  std::cout<<"The vector AX[i] is \n";
	  // show_vector(tempX);
	  tempaX=multiply(M, tempaX);
	  
	  AX(H[i][i], tempaX);
	  // std::cout<<"The vector lambda*MX[i] is:\n";
	  // show_vector(tempaX);
	  AYPX(-1.0, tempaX, tempX);
	  // std::cout<<"After AYPX() the vector AX[i] is :\n";
	  // show_vector(tempX);
	  // std::cout<<"the lambda*MX[i] is :\n";
	  // show_vector(tempaX);
	  // std::cout<<"The infi norm is : "<<infi_Norm(tempX)<<"\n";
	  if ( infi_Norm(tempX)>res)
	    {
	      res=infi_Norm(tempX);
	    }
	}
      
     
      if (res<tol)
	{
	  std::cout<<"After "<<k+1<<" iterations, \n";
	  lambda.clear();
	  for(int j=0; j<p;j++)
	    {
	      lambda.push_back(H[j][j]);
	      //lambda.push_back(1./R[j][j]); //for subspace method;
	    }
	  std::cout<<"The residual is "<<res<<"\n";
	  return;
	  //break;
	}
      k++;
    }
  lambda.clear();
  for(int j=0; j<p;j++)
    {
      //lambda.push_back(1./R[j][j]); //for subspace method;
      lambda.push_back(H[j][j]);
    }
  std::cout<<"The maximum number of iteration exceeded.\n";
  std::cout<<"The error is :: "<<res<<"\n";
  // std::cout<<"The numerical value of eigenvalues are :\n";
}

void EigenSolver::LOBPCGSolve(std::vector<std::vector<double>>& X,
			      std::vector<double>& lambda,
			      int p,
			      int max_iter,
			      double tol)
{
  std::vector<int> ind;
  std::vector<double> mu(p,0), tempx(A->m()), Ax, Mx, u(3*p,0), tempu;
  std::vector<std::vector<double>> res(A->m(),mu), tempM, tempMX, X0, Qk, tempU;
  double residual;
  int k=0;
  transpose(res);
  // std::cout<<"flag1\n";
  // std::cout<<"The matrix X is \n";
  // show_matrix(X);
  
  while (k < max_iter)
    {
      residual = 0;
      // std::cout<<"flag2\n";
      for(int i = 0;i < p; i++)
	{
	  for(int j=0;j<X.size();j++)
	    {
	      tempx[j]=X[j][i];
	    }
	  // std::cout<<"Temp x is \n";
	  // show_vector(tempx);
	  
	  Ax = multiply(A,tempx);
	  //  std::cout<<"Atempx is \n";
	  //  show_vector(Ax);
	  Mx = multiply(M,tempx);
	  // std::cout<<"Mtempx is \n";
	  //  show_vector(Mx);
	  //  std::cout<<"x'Mx is \n"<<inner(tempx,Mx)<<"\n";
	  //  std::cout<<"x'Ax is \n"<<inner(tempx,Ax)<<"\n";
	  mu[i] = inner(tempx,Mx)/inner(tempx,Ax);
	  AYPX(-1.0*mu[i], Ax, Mx);
	  res[i] = Mx;

	  if(infi_Norm(res[i])>residual)
	    {
	      residual = infi_Norm(res[i]);
	    }
	}
      
      // std::cout<<"The matrix res is \n";
      // show_matrix(res);
      //  std::cout<<"The vector mu is \n";
      //  show_vector(mu);
      //   std::cout<<"The residual is "<<residual <<"\n";
      //   std::cout<<"The tol is "<<tol<<"\n";
      //  std::cout<<"flag3\n";
      if(residual < tol)
	{
	  lambda = mu;
	  for(int i=0;i<lambda.size();i++)
	    {
	      lambda[i]=1.0/lambda[i];
	    }
	  std::cout<<"After "<<k+1<<" iterations, the methods converges\n";
	  return;
	}
      // std::cout<<"This is "<<k<<"-th iteration\n";
      if (k==0)
	{
	  tempM.clear();
	  //  std::cout<<"FLLLAG!\n";
	  for(int i=0;i<p;i++)
	    {
	      tempM.push_back(res[i]);
	    }
	  transpose(X);
	  for(int i=0;i<p;i++)
	    {
	      tempM.push_back(X[i]);
	    }
	  transpose(X);
	}
      else
	{
	  tempM.clear();
	  // std::cout<<"FFFFLAGGGG!!!\n";
	  for(int i=0;i<p;i++)
	    {
	      tempM.push_back(res[i]);
	    }
	  transpose(X);
	  for(int i=0;i<p;i++)
	    {
	      tempM.push_back(X[i]);
	    }
	  transpose(X);

	  transpose(X0);
	  for(int i=0;i<p;i++)
	    {
	      tempM.push_back(X0[i]);
	    }
	  transpose(X0);
	}
      transpose(tempM);
      
      // std::cout<<"flag4\n";
      // std::cout<<"The matrix tempM before GS is \n";
      // show_matrix (tempM);
      GS_A(tempM);
      //  std::cout<<"The matrix tempM is \n";
      //  show_matrix(tempM);
      tempMX = get_MX(tempM);
      // std::cout<<"The matrix tempMX is \n";
      //  show_matrix(tempMX);
      multiply(tempMX, tempM);
      //  std::cout<<"The matrix tempMX x tempM is \n";
      //  show_matrix(tempMX);
      
      identitymatrix(Qk, tempMX.size());
      
      QRSolver(tempMX, Qk, 1.e-6);
      // std::cout<<"The ritz value is \n";
      //  show_matrix(tempMX);
      //  std::cout<<"The ritz vector is \n";
      //  show_matrix(Qk);
      // multiply(X,Qk);

      //u.clear();
      for(int i=0;i<tempMX.size();i++)
	{
	  u[i] = tempMX[i][i];
	}
      ind.resize(u.size());

      iota(ind.begin(),ind.end(),0);

      sort(ind.begin(),ind.end(),[&u](int a,int b){return u[a]>u[b];});
      X0 = X;
      //multiply(X,Qk);
      tempU.clear();
      tempu.resize(Qk.size(),0);
      for(int i=0;i<p;i++)
	{
	  for(int j=0;j<Qk.size();j++)
	    {
	      tempu[j]=Qk[j][ind[i]];
	    }
	  tempU.push_back(tempu);
	}
      transpose(tempU);
      // std::cout<<"Matrix tempU is \n";
      // show_matrix(tempU);
      // std::cout<<"Matrix X is \n";
      multiply(tempM,tempU);
      X=tempM;
      //  show_matrix(X);
      k++;
    }

  std::cout<<"Maximum iteration numebr exceed.\n";
  lambda = mu;
  //std::vector<double> tmpzeros(lambda.size(),0);
  // To compute the smallest eigenvalues of A correponding to M 
  // which is equivalent to compute the largest eigenvalue of M
  // corresponding to A. Thus inverse operations is needed here.
  for(int i=0;i<lambda.size();i++)
    {
      lambda[i]=1.0/lambda[i];
    }
}
