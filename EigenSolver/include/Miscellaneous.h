#include<cmath>
#include<lac/sparsity_pattern.h>
#include<lac/sparse_matrix.h>

// Find the smallest index of the entry that absolute value is the infinity norm of the vector, i.e. |x_p| = || x || _ infinity.
void max_norm(std::vector<double> x, int &p)
{
  double max = 0;
  for (int i = 0;i < x.size(); i++)
    {
      if(fabs(x[i])>max)
	{
	  p = i;
	  max = fabs(x[i]);
	}
    }
}	

// Multiplication of a double parameter with a vector.
void AX(double a, std::vector<double>& x)
{
  for (int i = 0; i < x.size(); i++)
    {
      x[i]=a*x[i];
    }
}

// Multiplication of SparseMatrix A and vector x0;
std::vector<double> multiply(const dealii::SparseMatrix<double>* A, std::vector<double> x0)
{
  std::vector<double> x(x0.size(),0);
  //  std::cout<<"Flag 1!!\n";
  //std::cout<<A.n()<<"\n";
  
  if (A->n()!=x0.size())
    {
      std::cout<<"Function multiply() ERROR: The sizes of matrix and vector are not same! Please check it!\n";
    }
  else
    {
      // std::cout<<"TEST POINT 2 \n";
      for(int k=0;k<A->m();k++)
	{
	  //std::cout<<"this is the "<<k<<"th iterations \n";
		dealii::SparseMatrix<double>::const_iterator i=A->begin(k);
	  
	  while(i!=A->end(k))
	    {
	      x[k]+=i->value()*x0[i->column()];
	      ++i;
	    }
	}
      return x;
    }
}

// Compute x-ay and store the solution vector in x;
void AYPX(double a, std::vector<double> y, std::vector<double>& x)
{
  //std::vector<double> tempx(x.size());
  for(int i=0; i < x.size(); i++)
    {
      x[i] = a*y[i] + x[i];
    }

}

void show_vector(std::vector<double> a)
{
  for (int i=0;i<a.size();i++)
    {
      std::cout<<a[i]<<" ";
    }
  std::cout<<std::endl;
}

void show_matrix(std::vector<std::vector<double>> A)
{
  for(int i=0;i<A.size();i++)
    {
      for(int j=0;j<A[0].size();j++)
	{
	  std::cout<<A[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
      std::cout<<std::endl;
}
// Design for return transpose matrix of A;
void transpose(std::vector<std::vector<double>>& A)
{
  int n=A.size();
  int m=A[0].size();
    if (n < m)
      {
        for (int i = 0; i < n; i++)
	  {
            for (int j = 0; j < i; j++)
	      {
                double temp = A[j][i];
                A[j][i] = A[i][j];
                A[i][j] = temp;
	      }
	  }
        for (int i = n; i < m; i++)
	  {
            std::vector<double> temp;
            for (int j = 0; j < n; j++)
	      {
                temp.push_back(A[j][i]);
	      }
            A.push_back(temp);
	  }
        for (int i = 0; i < n; i++)
	  {
            for (int j = n; j < m; j++)
	      {
                A[i].pop_back();
	      }
	  }
      }
    else if(m<n)
      {
        for (int i = 0; i < m; i++)
	  {
            for (int j = 0; j < i; j++)
	      {
                double temp = A[j][i];
                A[j][i] = A[i][j];
		A[i][j] = temp;
	      }
	  }
        for (int i = 0; i < m; i++)
	  {
            for (int j = m; j < n; j++)
	      {
                A[i].push_back(A[j][i]);
	      }
	  }

        for (int i = m; i < n; i++)
	  {
            A.pop_back();
	  }

      }
    else
      {
	for(int i=0;i<n;i++)
	  {
	    for(int j=0;j<i;j++)
	      {
		double temp=A[j][i];
		A[j][i]=A[i][j];
		A[i][j]=temp;
	      }
	  }
      }

    //return A;

}


// this function will turn the matrix Q into the nxn identity matrix;
void identitymatrix(std::vector<std::vector<double>> &Q, int n)
{
  Q.clear();
  std::vector<double> temp(n,0);
  for (int i=0;i<n;i++)
    {
      Q.push_back(temp);
      Q[i][i]=1;
    }
}

//This function computes the matrix A times matrix B, and its result will be stored in A.
// It can be improved if reducing the use of the temp matrix to store the original values which
// is used in computing.
void multiply(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> B)
{
    int n = A.size();
    int m = B.size();
    int col = B[0].size();
    double tempX;
    std::vector<double> temp(col);
    std::vector<std::vector<double>> A0(n, temp);

    if(A[0].size()!=m)
      {
	std::cout<<"Function Multiply() ERROR: Please Check these two matrix size!!!\n";
	return;
      }
    /*   
    if (col < n)
      {
        for (int i = 0; i < n; i++)
	  {
            for (int j = 0; j < col; j++)
	      {
                double tempValue = 0;
                for (int t = 0; t < m; t++)
		  {
                    tempValue += A[i][t] * B[t][j];
		  }
                A0[i][j] = tempValue;
	      }
            for (int j = col; j < m; j++)
	      {
                A0[i].pop_back();
	      }
	  }

        A = A0;
      }
    else
      {
        for (int i = 0; i <n; i++)
	  {
            for (int j = 0; j < m; j++)
	      {
                double tempValue = 0;
                for (int t = 0; t < m; t++)
		  {
                    tempValue += A[i][t] * B[t][j];
		  }
                A0[i][j] = tempValue;
	      }
            for (int j = m; j < col; j++)
	      {
                double tempValue = 0;
                for (int t = 0; t < m; t++)
		  {
                    tempValue += A[i][t] * B[t][j];
		  }
                A0[i].push_back(tempValue);
	      }
	  }
        A = A0;
	}
    */
    for(int i=0;i<n;i++)
      {
	for(int j=0;j<col;j++)
	  {
	    tempX = 0;
	    for(int k=0;k<m;k++)
	      {
		tempX+=A[i][k]*B[k][j];
	      }
	    A0[i][j]=tempX;
	  }
      }
    A.clear();
    A=A0;
    //std::cout<<"The size of the matrix A is "<<A.size()<<" x "<<A[0].size()<<"\n";
}



// Computing matrix A times vector b, and using the A[0] to store the result and return it.
void multiply(std::vector<std::vector<double>> A, std::vector<double> &b)
{
  int m=A.size();
  int n=b.size();
  std::vector<double>tempx(m,0);
  for (int i=0;i<m;i++)
    {
      double temp=0;
      for(int j=0;j<n;j++)
	{
	  temp+=A[i][j]*b[j];
	}
      tempx[i]=temp;
    }
  b.clear();
  b=tempx;
  //return A[0];
}

// This function computes the infi-norm of vector x;
double infi_Norm(std::vector<double> x)
{
  double norm=0;
  for (int i=0;i<x.size();i++)
    {
      if (norm<fabs(x[i]))
	{
	  norm=fabs(x[i]);
	}
    }
  return norm;
}

//////////
// CG mehod;
double inner(std::vector<double> a, std::vector<double> b)
{
  int n=a.size();
  double sum=0;
  for(int i=0;i<n;i++)
    {
      sum+=a[i]*b[i];
    }
  return sum;
}

void CG(std::vector<std::vector<double>> A, std::vector<double>& x, std::vector<double> rhs, double tol=0.001, int max_iter=10)
{
  int n=A.size();
  std::vector<double> res(A.size(),0), Ap(A.size(),0), p(A.size(), 0);
  double delta=0, beta=0;
  //res=get_res(A,x,rhs); // res=b - Ax;
  for(int i=0;i<n;i++)
    {
      double temp=0;
      for(int j=0;j<n;j++)
	{
	  temp+=A[i][j]*x[j];
	}
      res[i]=rhs[i]-temp;
    }
  p=res;
  for(int i=0;i<n;i++)
    {
      p[i]=-res[i];
    }
  
  if (infi_Norm(p)<tol)
    {
      return;
    }

    //Ap=get_Ap(A,p);
  for(int i=0;i<n;i++)
    {
      double tempAp=0;
      for(int j=0;j<n;j++)
	{
	  tempAp+=A[i][j]*p[j];
	}
      Ap[i]=tempAp;
    }
  delta=inner(p,res)/inner(p,Ap);
  for (int i=0;i<n;i++)
    {
      x[i]=x[i]+delta*p[i];
    }
  int iter=0;
  
  // res=b - Ax;
  for(int i=0;i<n;i++)
    {
      double temp=0;
      for(int j=0;j<n;j++)
	{
	  temp+=A[i][j]*x[j];
	}
      res[i]=rhs[i]-temp;
    }
  while(infi_Norm(res)>tol&& iter < max_iter)
    {
      iter++;

      beta=inner(res, Ap)/inner(p, Ap);
      for(int k=0;k<p.size();k++)
        {
          p[k]=-res[k]+beta*p[k];
        }
      //get Ap;
      for(int i=0;i<n;i++)
	{
	  double tempAp=0;
	  for(int j=0;j<n;j++)
	    {
	      tempAp+=A[i][j]*p[j];
	    }
	  Ap[i]=tempAp;
	}
      delta=inner(p,res)/inner(p,Ap);

      //get new x;
      for (int i=0;i<n;i++)
	{
	  x[i]=x[i]+delta*p[i];
	}

      for(int i=0;i<n;i++)
	{
	  double temp=0;
	  for(int j=0;j<n;j++)
	    {
	      temp+=A[i][j]*x[j];
	    }
	  res[i]=rhs[i]-temp;
	}
    }
}

// Compute the Householder transports of matrix a and store the orthonormal matrix in X;
// Without any other requirement, the matrix X should be set as identity matrix.
void Householder(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> & X)
{
  //initialize the matrix X

  // Set tolerence to eliminate the situation that we get a number = - 0.0;
  double tolerence = 1.e-10;
  long double temp_1, temp_2;
  
  std::vector<double> temp(a.size(),0);
  // std::vector<std::vector<double>> tempX;
  // identity_matrix(tempX);
  //X.resize(a.size(),temp);// just for test this function, will be deleted after debug;
  // for(int j=0;j<X.size();j++)
  //  {
      //    X[j][j]=1;
  //  }
  int n=a.size();
  std::vector<std::vector<double>> Pk;//, tempX;
  // tempX = X;
  long double RSQ;
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
      if(q<tolerence)
        {
	  continue;
	}
      // Computing the value of alpha
      if(fabs(a[k+1][k])<tolerence)
        {
          alpha=-sqrt(q);
        }
      else
        {
          alpha=-sqrt(q)*a[k+1][k]/fabs(a[k+1][k]);
        }

      RSQ=alpha*alpha-alpha*a[k+1][k]; // RSQ = 2r^2;
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
          std::vector<double> tempttt(j-k,0);
          for(int i=k+1;i<=j;i++)
            {
              for(int t=k+1;t<n;t++)
                {
		  temp_1=v[i-k]*v[t-k]*Pk[t][j];
		  temp_2=(2*RSQ);
		  if (fabs(temp_1) < tolerence )
		    {
		      //temp_1=fabs(temp_1)
		      tempttt[i-k-1]+=0;
		      continue;
		    }
		  if(fabs(temp_2)<tolerence)
		    {
		      //temp_2=fabs(temp_2);
		      tempttt[i-k-1]+=0;
		      continue;
		    }
                  //tempttt[i-k-1]+=v[i-k]*v[t-k]*Pk[t][j]/(2*RSQ);
		  tempttt[i-k-1]+=1.e10*temp_1/(1.e10*temp_2);
                }
            }

          for (int i=k+1;i<=j;i++)
            {
              Pk[i][j]-=2*tempttt[i-k-1];
              Pk[j][i]=Pk[i][j];
            }
        }
      multiply(X,Pk);
     
    }
  // X=tempX;
}


//void TraceSolver::QR(std::vector<std::vector<double>> &a){};
///////////////
// This function computes the QR factorization of the tridiagonal matrix A.
// Input: tridiagonal matrix A, identity matrix Q
// Output: uptriangle matrix A_, the orthonormal matrix Q; A=Q*A_;
//
// tips: Maybe I can improve the function by replace the matrix A by two vectors an bn-1;
// an contains the diagonal entries of A; bn-1 contains the sub-diagonal entries of A;

void QR(std::vector<std::vector<double>> &a,   std::vector<std::vector<double>>& Q)
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



//Compute the shur decomposition of matrix a, and store the eigenvalues in the diagonal entries of a, store the eigenvectors in X.
// So if without any requirement, can set the initial X as a identity matrix.
void QRSolver(std::vector<std::vector<double>> &a, std::vector<std::vector<double>> &X, double tol=1.e-3)
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

      // std::cout<<"Before Householder a is :::;;\n";
      // show_matrix(a);
      // std::cout<<"X is :::::::;\n";
      // show_matrix(X);
      Householder(a, X);
      // std::cout<<"The matrix X after householder is :::::::::::::\n";
      // show_matrix(X);
      // std::cout<<"The matrix a after householder is ::::::::::::::\n";
      // show_matrix(a);

      identitymatrix(Qk,n);

      QR(a,Qk);

      // Compute the num-th iteration, get the Qk in this step;
      multiply(X,Qk);
      //             std::cout<<"The matrix X after QR is :::::::::::::\n";
      //             show_matrix(X);
      // compute R*Q beacuse I store the R(computed above) into A, So I directly use multiply
      // function to get A*Q into A=RQ.
      multiply(a,Qk);

      //update the entries of b, i.e. the subdiagonal entries of matrix A=RQ;
      for(int i=0;i<n-1;i++)
        {
          b.push_back(a[i+1][i]);
        }
    }
};


// This function using Gauss elimination to solve some simple linear equations;
void Gauss(std::vector<std::vector<double>> A, std::vector<double>&x, std::vector<double> b )
{
  double tol = 1.e-10, m;
  int n=A.size();
  std::vector<double> tempVec;
  for (int i=0;i<A.size();i++)
    {
      A[i].push_back(b[i]);
    }
  for(int i =0;i<n-1;i++)
    {
      int p = -1;
      // find p be the smallest integer with i<=p<=n and A[p,i]!=0;
      for(int k = i;k<n;k++)
        {
          if(fabs(A[k][i])>tol)
            {
              p = k;
              break;
            }
        }
      if (p==-1)
        {
          std::cout<<"No unique solution exists!\n";
          break;
        }

      if(p != i)
        {
          tempVec = A[p];
          A[p] = A[i];
          A[i] = tempVec;
        }
      for (int j=i+1;j<n;j++)
        {
          m=A[j][i]/A[i][i];
          AYPX(-1.0*m, A[i], A[j]);
        }
    }
  if (fabs(A[n-1][n-1])<tol)
    {
      std::cout<<"No unique solution exists!\n";
      return;
    }
  
  // Start backward substitution;
  x[n-1] = A[n-1][n]/A[n-1][n-1];
  for (int i=n-2; i >=0;i--)
    {
      double temp_sum = 0.;
      for(int j = i+1;j<n;j++)
	{
	  temp_sum+=A[i][j]*x[j];
	}
      x[i] = (A[i][n] - temp_sum)/A[i][i];
    }

}
