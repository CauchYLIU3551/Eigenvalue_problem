#include "elobpcg.h"
#include<lac/sparse_matrix.h>
#include<lac/sparsity_pattern.h>
elobpcg::elobpcg(){};
elobpcg::elobpcg(const Matrix& k, const Matrix& m)
{
  K=&k;
  M=&m;
  assert(K->m()==K->n());
  assert(M->m()==M->n());
}

elobpcg::elobpcg(const Matrix& k, const Matrix& m, const Matrix& e)
{
  K=&k;
  M=&m;
  E=&e;
  assert(K->m()==K->n());
  assert(M->m()==M->n());
  assert(E->m()==E->n());
}

std::vector<double> elobpcg::proj_K_shift(double shift, std::vector<std::vector<double>> Q, std::vector<double> u, std::vector<double> v)
{
  std::vector<std::vector<double>> Qt(Q);
  transpose(Qt);
  double delta=0;
  std::vector<double> Kv, Ku;
  Kv=multiply(K,v);
  Ku=multiply(K,u);
  std::vector<double> Qv, Qu;
  Qv=multiply(Qt,v);
  Qu=multiply(Qt,u);
  Qv=multiply(Q,Qv);
  Qu=multiply(Q,Qu);
  delta=(inner(u,Kv)+shift * inner(u,Qv))/(inner(u,Ku)+shift * inner(u,Qu)); // Compute the shift- K - inner product quotient to help compute the updated vector.
  for(int i=0;i<u.size();i++)
    {
      v[i]=delta*u[i];
    }
  return v;
}
  

void elobpcg::GS_K(double shift, std::vector<std::vector<double>> Q, std::vector<std::vector<double>>& X)
{
  transpose(X);
  std::vector<std::vector<double>> tempX(X);
  std::vector<std::vector<double>> Qt(Q);
  transpose(Qt);
  for (int i=1;i<tempX.size();i++)
    {
      for(int j=0;j<i;j++)
	{
	  std::vector<double> temp_proj;
	  temp_proj = proj_K_shift(shift, Q, X[j], tempX[i]); //compute the K-inner_product of X(:,j) and tempX(:,i);
	  double temp_proj_shift=0;
	  
	  for (int k=0;k<tempX[0].size();k++)
	    {
	      X[i][k] = X[i][k]-temp_proj[k];
	    }
	}
    }

  for (int k=0;k<X.size();k++)
    {
      std::vector<double> temp_Xk, temp_shift;
      temp_Xk = multiply(K,X[k]);
      temp_shift = multiply(Qt, X[k]);
      temp_shift = multiply(Q,temp_shift);

      double norm = 0;
      norm = inner(X[k], temp_Xk) + shift * inner(X[k], temp_shift);
      norm=sqrt(norm);
      for (int i=0; i<X[0].size();i++)
	{
	  X[k][i]=X[k][i]/norm;
	}
    }
  transpose(X);
}

std::vector<std::vector<double>> initial_Z(std::vector<std::vector<double>> Z)
{
  // X=Z(n+1:2n,:); Y=Z(1:n, :);
  int m=Z.size();
  int n=Z[0].size();
  assert(2*n==m);
  std::vector<std::vector<double>> X(n),Y(n);
  for (int i=0;i<n;i++)
    {
      X[i]=Z[i+n];
      Y[i]=Z[i];
    }
  std::vector<std::vector<double>> KX;
   std::vector<double> temp(n);
  // computing the matrix K*X=KX; and store it as (KX)* in fact;
  for (int j=0;j<n;j++)
    {
      for(int k=0;k<n;k++)
	{
	  temp[k]=X[k][j];
	}
      KX[j]=multiply(K,temp);
    }
  std::vector<std::vector<double>>XKX(KX);
  // computing (KX)* x X in fact = X* x (KX);
  multiply(XKX,X);
  // Computing cholesky factorization of XKX = L L', and compute L^-* to orthogonalize X;
  chol(XKX);
  low_solve(XKX);

  multiply(X, XKX);
  multiply(KX,XKX);

  // do similar commands on Y
  std::vector<std::vector<double>> MY;
  
  // computing the matrix M*Y=MY; and store it as (MY)* in fact;
  for (int j=0;j<n;j++)
    {
      for(int k=0;k<n;k++)
	{
	  temp[k]=Y[k][j];
	}
      MY[j]=multiply(M,temp);
    }
  std::vector<std::vector<double>>YMY(MY);
  // computing (KX)* x X in fact = X* x (KX);
  multiply(YMY,Y);
  // Computing cholesky factorization of XKX = L L', and compute L^-* to orthogonalize X;
  chol(YMY);
  low_solve(YMY);

  multiply(Y,  YMY);
  multiply(MY, YMY);

  std::vector<std::vector<double>> W(X);
  transpose(W);
  W=multiply(W,Y);
  // updating W as W=1/2 (W+ W');
  for (int i=0;i<W.size();i++)
    {
      for(int j=0;j<i;j++)
	{
	  W[i][j]=0.5*(W[i][j]+W[j][i]);
	  W[j][i]=W[i][j];
	}
    }

  // compute RX and RY;
}
