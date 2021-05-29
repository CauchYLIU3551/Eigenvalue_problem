#include<iostream>
#include<assert.h>
#include<lac/sparsity_pattern.h>
//#include<AFEPack/AMGSolver.h>
#include <trace/mintrace.h>
#include <trace/Miscellaneous.h>
/////////////////////////////////////////////
// These functions can be gathered into one .cpp file after testing;


/////////////////////////////////////////////
TraceSolver::TraceSolver(){};
TraceSolver::TraceSolver(const Matrix& a, const Matrix& m){
  A=&a;
  M=&m;
  assert(A->m()==A->n());
  assert(M->m()==M->n());
  assert(A->m()==M->m());

  // initializing the vector of theta;
  theta.resize(A->m());
}

void TraceSolver::rand_V(int p, std::vector<std::vector<double>>& V)
{
  //std::cout<<"Begin initialize V and X!\n";
  std::vector<double> tmp(p, 0);
  std::vector<std::vector<double>> tempV(A->m(), tmp);
  for (int i=0;i<p;i++)
    {
      //tempV[i][i]=sqrt(i+1)/(i+1);
      tempV[i][i]=1;
    }
  //V=tempV;
  //X=V;
  X=tempV;
  GS_M();
  // std::cout<<"This is matrix initial matrix X after GS_M()\n";
  //show_matrix(X);
  V=X;

};

// Computing the matrix Vt * A *V;
void TraceSolver::get_VtAV(std::vector<std::vector<double>>& V)
{
  int m=V.size();
  int n=V[0].size();
  std::vector<double> temp(n,0), tempAV(m,0);
  std::vector<std::vector<double>> tempV(n, temp);
  for(int i=0;i<n;i++)
    {
      for(int j=0;j<n;j++)
	{
	  double sum=0;
	  for(int k=0;k<m;k++)
	    {
	      tempAV[k]=V[k][j];
	    }
	  tempAV=multiply(A,tempAV);
	  for(int k=0;k<m;k++)
	    {
	      sum+=V[k][i]*tempAV[k];
	      tempV[i][j]=sum;
	    }
	}
    }
  V.clear();
  V=tempV;
};

// The function is used to compute the projection of v into u 
// under M-inner-productsComputing the project of <u, v> corresponding to M;
std::vector<double> TraceSolver::proj_M(std::vector<double> u, std::vector<double> v)
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
void TraceSolver::GS_M()
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

//void TraceSolver::Householder(std::vector<std::vector<double>> &a){};
void TraceSolver::Householder(std::vector<std::vector<double>> &a)
{
  //initialize the matrix X
  std::vector<double> temp(a.size(),0);
  // std::vector<std::vector<double>> tempX;
  // identity_matrix(tempX);
  //X.resize(a.size(),temp);// just for test this function, will be deleted after debug;
  // for(int j=0;j<X.size();j++)
  //  {
      //    X[j][j]=1;
  //  }
  int n=a.size();
  std::vector<std::vector<double>> Pk;

  for(int k=0;k<n-2;k++)
    {
      identitymatrix(Pk, n);
      double q=0;
      double alpha=0;
      // Computing the value of q in the Householder function;
      for(int i=k+1;i<n;i++)
	{
	  q+=a[i][k]*(a[i][k]);
	}

      // Computing the value of alpha
      if(a[k+1][k]==0)
	{
	  alpha=-sqrt(q);
	}
      else
	{
	  alpha=-sqrt(q)*a[k+1][k]/abs(a[k+1][k]);
	}
      
      double RSQ=alpha*alpha-alpha*a[k+1][k]; // RSQ = 2r^2;

      std::vector<double> v(n-k,0), u(n-k,0), z(n-k,0);
      v[0]=0;
      v[1]=a[k+1][k]-alpha;

      for(int j=2;j<v.size();j++)
	{
	  v[j]=a[k+j][k];
	}
      
      for(int j=0;j<n-k;j++)
	{
	  for(int i=1;i<n-k;i++)
	    {
	      u[j]+=a[k+j][k+i]*v[i];
	    }
	  u[j]/=RSQ;
	}
      
      double PROD=0;
      for(int i=1;i<n-k;i++)
	{
	  PROD+=v[i]*u[i];
	}
      
      for(int j=0;j<n-k;j++)
	{
	  z[j]=u[j]-PROD/(2*RSQ)*v[j];
	}

      for(int l=k+1;l<n-1;l++)
	{
	  for(int j=l+1;j<n;j++)
	    {
	      a[j][l]=a[j][l]-v[l-k]*z[j-k]-v[j-k]*z[l-k];
	      a[l][j]=a[j][l];
	    }
	  a[l][l]=a[l][l]-2*v[l-k]*z[l-k];
	}
      
      a[n-1][n-1]=a[n-1][n-1]-2*v[n-k-1]*z[n-k-1];
	
      for(int j=k+2;j<n;j++)
	{
	  a[k][j]=0;
	  a[j][k]=0;
	}

      a[k+1][k]=a[k+1][k]-v[1]*z[0];
      a[k][k+1]=a[k+1][k];


      for (int j=k+1;j<n;j++)
	{
	  std::vector<double> temp(j-k,0);
	  for(int i=k+1;i<=j;i++)
	    {
	      for(int t=k+1;t<n;t++)
		{
		  temp[i-k-1]+=v[i-k]*v[t-k]*Pk[t][j]/(2*RSQ);
		}
	    }
	  for (int i=k+1;i<=j;i++)
	    {
	      Pk[i][j]-=2*temp[i-k-1];
	      Pk[j][i]=Pk[i][j];
	    }
	}
      multiply(X,Pk);
    }

}

//void TraceSolver::QR(std::vector<std::vector<double>> &a){};
///////////////
// This function computes the QR factorization of the tridiagonal matrix A.
// Input: tridiagonal matrix A, identity matrix Q
// Output: uptriangle matrix A_, the orthonormal matrix Q; A=Q*A_;
//
// tips: Maybe I can improve the function by replace the matrix A by two vectors an bn-1;
// an contains the diagonal entries of A; bn-1 contains the sub-diagonal entries of A;

void TraceSolver::QR(std::vector<std::vector<double>> &a,   std::vector<std::vector<double>>& Q)
{
  int n=a.size();
  // std::cout<<"The matrix Qk before QR is ::::::::::!!!!!!!!!!!\n";
  // show_matrix(Q);
  //   std::cout<<"The matrix a before QR is ::::::::::!!!!!!!!!!!\n";
  //  show_matrix(a);

  for (int i=0;i<a.size()-1;i++)
    {
      double temp_theta=0;
      double tempa1=a[i][i],tempa2=a[i+1][i+1],tempb=a[i+1][i],tempc=a[i][i+1];
      temp_theta=atan(a[i+1][i]/a[i][i]);
      double c=cos(temp_theta), s=sin(temp_theta);

      a[i][i]=c*tempa1+s*tempb;
      a[i][i+1]=c*tempc+s*tempa2;
      a[i+1][i+1]=-s*tempc+c*tempa2;
      a[i+1][i]=-s*tempa1+c*tempb;
      if (i<a.size()-2)
	{
	  a[i][i+2]=s*a[i+1][i+2];
	  a[i+1][i+2]=c*a[i+1][i+2];
	}


      // update the orthogonal matrix Q by Q_k=Q_k-1*Q;
      
      for(int k=0;k<n;k++)
	{
	  double tmp1=Q[k][i],tmp2=Q[k][i+1];
	  Q[k][i]=Q[k][i]*c+Q[k][i+1]*s;
	  Q[k][i+1]=tmp1*(-s)+tmp2*c;
	}
      
    }
  // std::cout<<"The matrix Qk in QR is ::::::::::!!!!!!!!!!!\n";
  // show_matrix(Q);
  //   std::cout<<"The matrix aQk in QR is ::::::::::!!!!!!!!!!!\n";
  // show_matrix(a);
  //multiply(X,Q); this update the X=V*Qk but in QR factorization, it does not contain this command;
  //multiply(a,Q); this update a=RQ, because after this QR function, a=R, after this command, a=RQ;
}

void TraceSolver::QRSolver(std::vector<std::vector<double>> &a, double tol)
{
  int n=a.size();
  // Qk as the initial matrix for every QR factorization, multiplying it together to get the
  // final unitary matrix;
  std::vector<std::vector<double>> Qk;
  
  // Get the subdiagonal entries of the matrix A and verify if its norm smaller than the tolerance;
  // If it is small enough that means we diagonalize the matrix A and the diagonal entries are
  // eigenvalues of A.
  std::vector<double> b(n-1,1);
  for(int i=0;i<n-1;i++)
    {
      b.push_back(a[i+1][i]);
    }

  int num=0;
  while (infi_Norm(b)>tol&&num<100)
    {
      num++;

      b.clear();

      Householder(a);
      // std::cout<<"The matrix X after householder is :::::::::::::\n";
      //    show_matrix(X);
      
      identitymatrix(Qk,n);

      QR(a,Qk);
     
      // Compute the num-th iteration, get the Qk in this step;
      multiply(X,Qk);
      //	     std::cout<<"The matrix X after QR is :::::::::::::\n";
      //	     show_matrix(X);
      // compute R*Q beacuse I store the R(computed above) into A, So I directly use multiply
      // function to get A*Q into A=RQ.
      multiply(a,Qk);

      //update the entries of b, i.e. the subdiagonal entries of matrix A=RQ;
      for(int i=0;i<n-1;i++)
	{
	  b.push_back(a[i+1][i]);
	}
    }
  theta.clear();
  theta.resize(a.size());
  for(int i=0; i<a.size(); i++)
    {
      theta[i]=a[i][i];
    }
  /*
  std::cout<<"This is the matrix a\n";
  show_matrix(a);*/
};

void TraceSolver::get_MX()
{
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
};


int globaliter=0, global_i=0;
void TraceSolver::get_Px(std::vector<double> & x)
{
  //std::cout<<"Entering get_Px()!!!!!!!!!!!!!!!!!!!!!!\n";
  std::vector<double>temp_x(x), MXx;
  // std::cout<<"This x #!##!#!#!#!#!#!#\n";
  //  show_vector(x);
    //std::cout<<"This is MX\n";
  /*
  if(globaliter==720&&global_i==1)
    {
      std::cout<<"The X is !!!!!!!!!!!!!!!\n";
      show_matrix(X);
      //get_MX();

      std::vector<double> tempx2(M->n()),TempMx;
      for(int k=0;k<M->m();k++)
	{
	  tempx2[k]=1;
	  TempMx=multiply(M,tempx2);
	  for(int q=0;q<TempMx.size();q++)
	    {
	      std::cout<<TempMx[q]<<" ";
	    }
	  std::cout<<std::endl;
	  tempx2[k]=0;
	}
      std::cout<<"The MX is \n";
      show_matrix(MX);
      std::cout<<"The vector AXei is ::\n";
      for(int i=0;i<x.size();i++)
	{
	  std::cout<<x[i]<<" ";
	}
      std::cout<<std::endl;
    }
  */
  
      
  multiply(MX,temp_x);

  /*
  if(globaliter==720&&global_i==1)
    {
      std::cout<<"The vector (MX)' * AXei is :::::::::::::\n";
      for(int i=0;i<tempx.size();i++)
	{
	  std::cout<<tempx[i]<<" ";
	}
      std::cout<<std::endl;
      
    }
    //  std::cout<<"This is tempx times MX\n";
    //show_vector(tempx);
// computing the (MX)t * (MX);
  //  std::cout<<"This is matrix MX\n";
  //  show_matrix(MX);
  */
  std::vector<double> temp(MX.size());
  std::vector<std::vector<double>> MXtMX(MX.size(),temp);
  for (int i=0;i<MX.size();i++)
    {
      for(int j=0;j<MX.size();j++)
	{
	  double sum=0;
	  for(int k=0;k<MX[0].size();k++)
	    {
	      sum+=MX[i][k]*MX[j][k];
	    }
	  MXtMX[i][j]=sum;
      	}
    }
  //  std::cout<<"This is MXtMX\n";
  //show_matrix(MXtMX);
  /*if(globaliter==720&&global_i==1)
    {
      std::cout<<"This is MXtMX\n";
      show_matrix(MXtMX);
      std::cout<<"\n";
    }*/

  // solve the equation: (MXtMX)^-1 x = b, i.e. (MXtMX) b = x to get vector b;
  std::vector<double> RHS(temp_x);
  temp_x.clear();
  temp_x.resize(RHS.size(),0);
  // std::cout<<"This is tempx.resize()\n";
  //show_vector(tempx);
  // std::cout<<"This is vector RHS\n";
  //show_vector(RHS);
  double tol2=1.0e-5;
  CG(MXtMX, temp_x, RHS, tol2, temp_x.size());
  for (int i=0;i<x.size();i++)
    {
      double sum=0;
      for(int j=0;j<MX.size();j++)
	{
	  sum+=MX[j][i]*temp_x[j];
	}
      x[i]=x[i]-sum;
    }

  // std::cout<<"This is vector x after compute!!!\n";
  // show_vector(x);
  //std::cout<<std::endl;
    //<<std::endl;
}

void TraceSolver::get_Ap(std::vector<double> p)
{
  // A great idea: while compute Ap, exactly (PAP)p here. There is a inverse
  // matrix in P, but it is unnecessary to compute the inverse matrix!
  // For B^-1 * x = b, while we know B and x, we can solve the equation:
  // B* b = x, to get the vector b!
  //std::cout<<"This is TraceSolver::get_Ap()\n";
  
  get_Px(p);
  //std::cout<<"This is Px\n";
  //show_vector(p);
  p=multiply(A, p);
  get_Px(p);
  Ap.reinit(p.size());
  for(int i=0;i<Ap.size();i++)
    {
      Ap[i]=p[i];
    }
};

void TraceSolver::get_res(std::vector<double> x, std::vector<double> r)
{
  //std::cout<<"This is TraceSolver::get_res()\n";

  res.reinit(x.size());
  //std::cout<<"This is x\n";
  // show_vector(x);
  std::vector<double> tempx(x);
  //get_Px(tempx);
  get_Px(x);
  // std::cout<<"This is Px\n";
  //show_vector(tempx);
  // show_vector(x);
  tempx=multiply(A,x);
  x=tempx;
  // std::cout<<"This is APx\n";
  // show_vector(x);
  get_Px(x);
  // std::cout<<"This is PAPx\n";
  //  show_vector(x);
  for(int i=0;i<res.size();i++)
    {
      res[i]=r[i]-x[i];
    }

  //std::cout<<"Finishing TraceSolver::get_res()\n";
}

// Computing the max error between AX and MXtheta;
double TraceSolver::get_residual()
{
  double residual=0,tempres;
  std::vector<std::vector<double>> tempX(X);
  transpose(tempX);
  std::vector<double> tempAX,tempMXtheta;
  
  for(int i=0;i<tempX.size();i++)
    {
      tempAX=multiply(A,tempX[i]);
      tempMXtheta=MX[i];
      for(int j=0;j<MX[i].size();j++)
	{
	  //tempMXtheta[j]=tempMXtheta[j]*theta;
	  if(residual<fabs(tempAX[j]-tempMXtheta[j]*theta[i]))
	    {
	      residual=fabs(tempAX[j]-tempMXtheta[j]*theta[i]);
	    }
	}
      //tempMXtheta=multiply(M,)
      //get_MX();
    }
  return residual;
}

/*void TraceSolver::solve(std::vector<double>& x,double tol, u_int step){}*/

// solve function computes the smallest p eigenvalues of the Matrix A and M;
// and the eigenvectors will be stored in the matrix X;
void TraceSolver::mintrace(int p, double tol, u_int max_iter)
{
  std::vector<std::vector<double>> V(A->n());
  std::vector<std::vector<double>> tempX;

  // using rand_V function to get the initial matrix V which is orthogonal corresponding with M;
  rand_V(p, V);// at this time it will send X = V; cause X=VU;
  //get_MX();
  int iter = 0;
  while(iter<max_iter)
    {
      // get_MX();
      // std::cout<<"The matrix X is :::::::::::::\n";
      //  show_matrix(X);

      /*
      if(iter==720)
	{
	  //identitymatrix(X,A->m());
	  std::cout<<"The matrix V before get_VtAV is :::::::::::::\n";
	  show_matrix(V);
	  }*/
      get_VtAV(V);
      /*
      if(iter==720)
	{
	  std::cout<<"The matrix VAV is :::::::::::::\n";
	  show_matrix(V);
	  }*/
      
      QRSolver(V);// return the diagonal entries of VtAV in theta and the V_k* U_k stored in X;
      /*
      if(iter==720)
	{
	  std::cout<<"The matrix VAV after QRSolver is :::::::::::::\n";
	  show_matrix(V);
	  std::cout<<"\n Meanwhile, the X is :::::::\n";
	  show_matrix(X);
	  }*/
      //Householder(V);//this function will return the tridiagonal function in V and the Householder matrices stored in X;
      //QR(V);
      //   std::cout<<"The matrix X is :::::::::::::\n";
      //    show_matrix(X);
      tempX=X;
      get_MX();

      /*
      if(iter==720)
	{
	  std::cout<<"The matrix MX is :::::::::::::\n";
	  show_matrix(MX);
	  }*/
	    
      double residual=get_residual();
      if (residual<tol)
	{
	  //std::cout<<"The iter number is : "<<iter<<"\n";
	  // std::cout<<"At this time MX is :\n";
	  // show_matrix(MX);
	  //std::cout<<"\n";
	  //std::cout<<"The theta is :\n";
	  //for(int i=0;i<theta.size();i++)
	  //{
	  //  std::cout<<theta[i]<<" ";
	  // }
	  /// std::cout<<"\n";
	  std::cout<<"The distinction error is :\n"
		   <<residual<<std::endl;
	  break;
	}

      /*
      if(iter==720)
	{
	  std::cout<<"The matrix X before CG solution process is :::::::::::::\n";
	  show_matrix(X);
	  }*/

      ///////////////////
      /////  globaliter=iter;
      //////////////////
      
      for (int i=0;i<p;i++)
	{
	  // std::cout<<"The matrix X is :::::::::::::\n";
	  //  show_matrix(X);

	  //
	  ////// global_i=i+1;
	  //
	  
	  std::vector<double> delta(A->m(),0), rhs(A->m(),0);// using a n x 1 dimension initial vector to compute the Conjugate gradient process.
	  //computing rhs vector!!
	  for(int j=0;j<rhs.size();j++)
	    {
	      rhs[j]=X[j][i];
	    }
	  /*
	  if(iter==720)
	    {
	      std::cout<<"X*e::::::::::"<<i+1<<std::endl;
	      //   show_vector(rhs);
	      for(int i=0;i<rhs.size();i++)
		{
		  std::cout<<rhs[i]<<" ";
		}
	      std::cout<<std::endl;
	    }*/
	  rhs=multiply(A, rhs);
	  //std::cout<<"This is A*X[i]\n";

	  /*
	  if(iter==720)
	    {
	      std::cout<<"AX*e::::::::::"<<i+1<<std::endl;
	      //   show_vector(rhs);
	      for(int i=0;i<rhs.size();i++)
		{
		  std::cout<<rhs[i]<<" ";
		}
	      std::cout<<std::endl;
	      }*/
		  
	     //  show_vector(rhs);
	  get_Px(rhs);

	  ///
	  /////////////global_i=0;
	  ///
	  
	  //std::cout<<"This is P*A*X[i]\n";

	  /*
	  	  if(iter==720)
	    {
	      std::cout<<"PAX*e::::::::::"<<i+1<<std::endl;
	      //   show_vector(rhs);
	      for(int i=0;i<rhs.size();i++)
		{
		  std::cout<<rhs[i]<<" ";
		}
	      std::cout<<std::endl;
	      }*/
		  
	  // show_vector(rhs);

		  /*	  
	  if(iter==720)
	    {
	      std::cout<<"This is the initial vector delta for "<<i+1<<" test!!@#@!#@!#@!#@!#@!#!@#\n";
	      show_vector(delta);
	      }*/
	  solve(delta,rhs,1.0e-2,A->m());
	  // std::cout<<"flag1\n";

	  /*
	  if(iter == 720)	  
	    {
	      std::cout<<"This is the solution vector delta\n";
	      show_vector(delta);
	      std::cout<<"\n !!!\n";
	      }*/
	  // using delta and X to compute V_k+1 stored in V and X, then do iteration
	  // again;
	  // G-S orthogonormalize function to make V_k+1 is ortho by M;
	  //
	  // update the i th column of X;
	  for(int k=0;k<delta.size();k++)
	    {
	      X[k][i]=X[k][i]-delta[k];
	    }
	  // std::cout<<"This is the "<<i<<" th iteration\n";
	}
     
      ///  std::cout<<"This is the X_k+1\n";
      ////  show_matrix(X);


      /*
      if(iter==720)
	{
	  std::cout<<"The matrix X-delta is :::::::::::::\n";
	  show_matrix(X);
	  }
      */
      GS_M();
    /*
      if(iter==720)
	{
	  std::cout<<"The matrix X after GS_M is :::::::::::::\n";
	  show_matrix(X);
	  }*/
      /////std::cout<<"This is the matrix X after GS_M\n";
      // show_matrix(X);
      V=X;
      //std::cout<<"The "<<iter+1<<" th iteration finished!!!!!!!!!!!!!!!!!\n";
      iter++;
      //if(iter==720)
      //	{
      //	  std::cout<<"Print the matrix V_k\n";
      //  show_matrix(V);
      //	}
      
    }
  // std::cout<<"The final matrix X is :\n";
  // show_matrix(tempX);
  // std::cout<<"The final matrix V is :\n";
  // show_matrix(V);
  std::cout<<"After "<<iter<<" th iteration, the "<<p<<" smallest eigenvalues of the equation are stored in lambda member.\n";
  lambda.clear();
  lambda.resize(p);
  transpose(tempX);
  //tempX=X;
  for(int i=0;i<p;i++)
    {
      std::vector<double> temp_AX;
      temp_AX=multiply(A,tempX[i]);
      lambda[i]=inner(tempX[i],temp_AX);
    }
  
}



