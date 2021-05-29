#include<iostream>
#include<fstream>
#include<vector>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>

using namespace dealii;

// Computing matrix A times vector b, and using the A[0] to store the result and return it.
std::vector<double> multiply(std::vector<std::vector<double>> A, std::vector<double> b)
{
  int m=A.size();
  int n=b.size();
  for (int i=0;i<m;i++)
    {
      double temp=0;
      for(int j=0;j<n;j++)
	{
	  temp+=A[i][j]*b[j];
	}
      A[0][i]=temp;
      //std::cout<<temp;
      //std::cout<<"The result of "<<i+1<<"th iteration\n";			
      
    }

  return A[0];
}

// Computing u'*v;
double multiply(std::vector<double> u, std::vector<double> v)
{
  if(u.size()!=v.size())
    {
      std::cout<<"Multiply() ERROR: u'*v error! Please Check the size of u and v!";
    }
  else
    {
      double sum=0;
      for(int i=0;i<u.size();i++)
	{
	  sum+=u[i]*v[i];
	}
      
      return sum;
    }
}

// Design for return transpose matrix of A;
std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> A)
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

    return A;

}

//This function computes the matrix A times matrix B, and its result will be stored in A.
// It can be improved if reducing the use of the temp matrix to store the original values which
// is used in computing.
void multiply(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> B)
{
    int n = A.size();
    int m = B.size();
    int col = B[0].size();
    std::vector<double> temp(m);
    std::vector<std::vector<double>> A0(n, temp);

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

// This function computes the infi-norm of vector x;
double infi_Norm(std::vector<double> x)
{
  double norm=0;
  for (int i=0;i<x.size();i++)
    {
      if (norm<abs(x[i]))
	{
	  norm=abs(x[i]);
	}
    }
  return norm;
}

//////
// Multiplication of SparseMatrix A and vector x0;
std::vector<double> multiply(const SparseMatrix<double> &A, std::vector<double> x0)
{
  std::vector<double> x(x0.size(),0);
  //  std::cout<<"Flag 1!!\n";
  //std::cout<<A.n()<<"\n";
  
  if (A.n()!=x0.size())
    {
      std::cout<<"Function multiply() ERROR: The sizes of matrix and vector are not same! Please check it!\n";
    }
  else
    {
      // std::cout<<"TEST POINT 2 \n";
      for(int k=0;k<A.m();k++)
	{
	  //std::cout<<"this is the "<<k<<"th iterations \n";
	  SparseMatrix<double>::const_iterator i=A.begin(k);
	  
	  while(i!=A.end(k))
	    {
	      x[k]+=i->value()*x0[i->column()];
	      ++i;
	    }
	}
      return x;
    }
}


//void Householder(std::vector<std::vector<double>> *A, std::vector<std::vector<double>> *P)
//void Householder(std::vector<std::vector<double>> *A, std::vector<std::vector<double>> *P)

// This function computes the householder process of the symmetric definite matrix, it return
// two matrices. The tridiagonal matrix coordinating with A and the hermitian unitary matrix P.
void Householder(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &P)
{
  int n=A.size();
  std::vector<std::vector<double>> Pk;

  for(int k=0;k<n-2;k++)
    {
      identitymatrix(Pk, n);
      double q=0;
      double alpha=0;
      // Computing the value of q in the Householder function;
      for(int i=k+1;i<n;i++)
	{
	  q+=A[i][k]*(A[i][k]);
	}

      // Computing the value of alpha
      if(A[k+1][k]==0)
	{
	  alpha=-sqrt(q);
	}
      else
	{
	  alpha=-sqrt(q)*A[k+1][k]/abs(A[k+1][k]);
	  // std::cout<<"the abs A[k+1][k] is "<<A[k+1][k]<<std::endl;
	}
      
      double RSQ=alpha*alpha-alpha*A[k+1][k]; // RSQ = 2r^2;

      //std::cout<<"Check Point 2, This is RSQ:::"<<RSQ<<"\n";
	
      std::vector<double> v(n-k,0), u(n-k,0), z(n-k,0);
      v[0]=0;
      v[1]=A[k+1][k]-alpha;

      //std::cout<<"Check Point 3::: "<<v[1]<<" \n";
      // w=(1/sqrt(2*RSQ)*v=1/2r*v);
      for(int j=2;j<v.size();j++)
	{
	  v[j]=A[k+j][k];
	}
      
      // std::cout<<"CheckPPPPPoint555 \n";
      
      for(int j=0;j<n-k;j++)
	{
	  for(int i=1;i<n-k;i++)
	    {
	      u[j]+=A[k+j][k+i]*v[i];
	    }
	  u[j]/=RSQ;
	}

      //std::cout<<"Check POOOIONIOHIHO\n";
      
      double PROD=0;
      for(int i=1;i<n-k;i++)
	{
	  PROD+=v[i]*u[i];
	}

      //std::cout<<"CheckPoint 4!!!!!!!!!\n";
      
      for(int j=0;j<n-k;j++)
	{
	  z[j]=u[j]-PROD/(2*RSQ)*v[j];
	}

      for(int l=k+1;l<n-1;l++)
	{
	  for(int j=l+1;j<n;j++)
	    {
	      A[j][l]=A[j][l]-v[l-k]*z[j-k]-v[j-k]*z[l-k];
	      A[l][j]=A[j][l];
	    }
	  A[l][l]=A[l][l]-2*v[l-k]*z[l-k];
	}
      
      A[n-1][n-1]=A[n-1][n-1]-2*v[n-k-1]*z[n-k-1];
	
      for(int j=k+2;j<n;j++)
	{
	  A[k][j]=0;
	  A[j][k]=0;
	}

      A[k+1][k]=A[k+1][k]-v[1]*z[0];
      A[k][k+1]=A[k+1][k];


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
      
	  //temp=multiply(W,Pj);// to get the jth column of the W*P;
	  for (int i=k+1;i<=j;i++)
	    {
	      Pk[i][j]-=2*temp[i-k-1];
	      Pk[j][i]=Pk[i][j];
	    }
	}
      multiply(P,Pk);

    }

}


///////
// This function computes the QR factorization of the tridiagonal matrix A.
// Input: tridiagonal matrix A, identity matrix Q
// Output: uptriangle matrix A_, the orthonormal matrix Q; A=Q*A_;
//
// tips: Maybe I can improve the function by replace the matrix A by two vectors an bn-1;
// an contains the diagonal entries of A; bn-1 contains the sub-diagonal entries of A;

void QR(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &Q)
{
  int n=A.size();

  // obtain the diagonal and subdiagonal entries of matrix A;
  /*
  std::vector<double>an(n),bn(n-1);
  for (int i=0;i<A.size()-1;i++)
    {
      an[i]=A[i][i];
      bn[i]=A[i+1][i];
    }
  an[n-1]=A[n-1][n-1];
  */
  
  for (int i=0;i<A.size()-1;i++)
    {
      double theta=0;
      double tempa1=A[i][i],tempa2=A[i+1][i+1],tempb=A[i+1][i],tempc=A[i][i+1];
      theta=atan(A[i+1][i]/A[i][i]);
      double c=cos(theta), s=sin(theta);

      A[i][i]=c*tempa1+s*tempb;
      A[i][i+1]=c*tempc+s*tempa2;
      A[i+1][i+1]=-s*tempc+c*tempa2;
      A[i+1][i]=-s*tempa1+c*tempb;
      if (i<A.size()-2)
	{
	  A[i][i+2]=s*A[i+1][i+2];
	  A[i+1][i+2]=c*A[i+1][i+2];
	}


      // update the orthogonal matrix Q by Q_k=Q_k-1*Q;
      
      for(int k=0;k<n;k++)
	{
	  double tmp1=Q[k][i],tmp2=Q[k][i+1];
	  Q[k][i]=Q[k][i]*c+Q[k][i+1]*s;
	  Q[k][i+1]=tmp1*(-s)+tmp2*c;
	}
      
    }

}




// This function computes the eigenvalue of the matrix A by QR methods, and it will also return the
// matrix Q during the process;
void QRSolver(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &Q, double tol=0.001)
{
  int n=A.size();
  // using householder transform matrix A into the tridiagonal matrix
  // store it in A and store the Householder unitary matrix in Q; A=QA_Q; 
  //Householder(A,Q);

  // Qk as the initial matrix for every QR factorization, multiplying it together to get the
  // final unitary matrix;
  std::vector<std::vector<double>> Qk;
  
  // Get the subdiagonal entries of the matrix A and verify if its norm smaller than the tolerance;
  // If it is small enough that means we diagonalize the matrix A and the diagonal entries are
  // eigenvalues of A.
  std::vector<double> b(n-1,1);
  for(int i=0;i<n-1;i++)
    {
      b.push_back(A[i+1][i]);
    }

  int num=0;
  while (infi_Norm(b)>tol&&num<1000)
    {
      num++;

      b.clear();

      Householder(A,Q);
      
      identitymatrix(Qk,n);

      QR(A,Qk);
     
      // Compute the num-th iteration, get the Qk in this step;
      multiply(Q,Qk);
      // compute R*Q beacuse I store the R(computed above) into A, So I directly use multiply
      // function to get A*Q into A=RQ.
      multiply(A,Qk);

      //update the entries of b, i.e. the subdiagonal entries of matrix A=RQ;
        for(int i=0;i<n-1;i++)
	  {
	    b.push_back(A[i+1][i]);
	  }
    }
}


class MinTrace
{
 public:
  MinTrace();
  MinTrace(SparseMatrix<double> A, SparseMatrix<double> M, int p=1);

  // get the random initial matrix V which satisifies that: V'*M*V=I_p;
  std::vector<std::vector<double>>rand_V(); 

  // compute the orthonormal vectors corresponding to M from the original vectors U={u_1,u_2,...,u_p};
  std::vector<std::vector<double>> M_GS(std::vector<std::vector<double>> U);
  
  // Compute the solution vector of the equation of X_k[i] and PAP*d_i=PA*X_k[i];
  std::vector<double>CG(std::vector<double> X_ki);
  // get the p smallest eigenvalues of A correponding with M;
  std::vector<double> min_trace(int max_iter, double eps);
  
  
 private:
  SparseMatrix<double> A;
  SparseMatrix<double> M;
  int p; // The number of the eigenvalues to be solved;
};

MinTrace::MinTrace()
{
  int p=3;
  SparsityPattern sparsity_pattern(10,10,{2,3,3,3,3,3,3,3,3,2});
  SparsityPattern sp2(10,10,1);
  //sparsity_pattern.add(1,2);
  sparsity_pattern.add(0,1);
  sparsity_pattern.add(9,8);
  for (int i=1;i<9;i++)
    {
      sparsity_pattern.add(i,i+1);
      sparsity_pattern.add(i,i-1);
    }
  
  A.reinit(sparsity_pattern);
  M.reinit(sp2);

  for (int k=0;k<A.m();k++)
    {
      SparseMatrix<double>::iterator i=A.begin(k);
      i->value()=2;
      while(i!=A.end(k))
	{
	  i->value()+=1;
	  ++i;
	}
    }

  SparseMatrix<double>::iterator i=M.begin();
  double num=0;
  while(i!=M.end())
    {
      num=num+1;
      i->value()=num;
      ++i;
    }
}

// ATTENTION!!! There is a difficulty that SparseMatrix can not be copied directly!
MinTrace::MinTrace(SparseMatrix<double> A0, SparseMatrix<double> M0, int p0)
{
  
};

std::vector<std::vector<double>> MinTrace::M_GS(std::vector<std::vector<double>> U)
{
  std::vector<std::vector<double>> V;
  V=transpose(U);
}

std::vector<std::vector<double>> MinTrace::rand_V()
{
  std::cout<<"CheckPoint1\n";
  std::vector<double> x(M.n(),0);
  for (int i=0;i<x.size();i++)
    {
      std::cout<<x[i]<<" ";
    }
  std::cout<<"\n";
  std::cout<<"CheckPoint666 \n";
  std::vector<std::vector<double>> V(p,x);
  std::cout<<"CheckPoint12314113\n";
  
  V[0][0]=1;
  V[1][1]=sqrt(2.0)/2;
  V[2][2]=sqrt(3.0)/3;
  std::cout<<"231wqrdasd24121321321321adsadq2q1312321sada\n";
  return V;
}

std::vector<double> MinTrace::min_trace(int max_iter=10, double eps=1e-03)
{
  std::cout<<M.n()<<"   it is a test!\n";
  
  std::vector<double> x(p); // define the eigenvector;
  std::vector<std::vector<double>> V; // V is an n x p matrix and stored as V[p][n];

  V=rand_V();
  std::cout<<"wewqewqewqewqeqw\n";
  //  for(int k=0;k<max_iter;k++)
    for(int k=0;k<1;k++)
    {
      std::cout<<"CheckPoint222\n";
      std::vector<std::vector<double>> VMV(p,x),tmpVMV(p); // VMV is a p x p matrix;
      std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAA\n";
      for(int i=0;i<p;i++)
	{
	  tmpVMV[i]=multiply(M,V[i]);
	  std::cout<<"CheckPoint111111111\n";
	}

      std::cout<<"CheckPoint2\n";
      for(int i=0;i<p;i++)
	{
	  for(int j=0;j<p;j++)
	    {
	      double tmp=0;
	      for(int l=0;l<M.n();l++)
		{
		  tmp+=V[i][l]*tmpVMV[j][l];
		}
	      VMV[i][j]=tmp;
	    }
	}

      for(int i=0;i<p;i++)
	{
	  for(int j=0;j<p;j++)
	    {
	      std::cout<<VMV[i][j]<<" ";
	    }
	  std::cout<<"\n";
	}
    }
  
   return x;
}

int main()
{
  SparseMatrix<double> A, M;
  // std::ofstream sparsematrix1 ("original_matrix.1");
  //A.print(sparsematrix1);
  std::vector<unsigned int> row_length(10,10);
  unsigned int t=3;
  SparsityPattern sparsity_pattern(10,10,{2,3,3,3,3,3,3,3,3,2});
  SparsityPattern sp2(10,10,1);
  //sparsity_pattern.add(1,2);
  sparsity_pattern.add(0,1);
  sparsity_pattern.add(9,8);
  for (int i=1;i<9;i++)
    {
      sparsity_pattern.add(i,i+1);
      sparsity_pattern.add(i,i-1);
    }
  
  A.reinit(sparsity_pattern);
  M.reinit(sp2);

  // In this way, I can construct a SparseMatrix to be the test data for the algorithm.
  for (int k=0;k<A.m();k++)
    {
      SparseMatrix<double>::iterator i=A.begin(k);
      i->value()=2;
      while(i!=A.end(k))
	{
	  i->value()+=1;
	  ++i;
	}
    }

  SparseMatrix<double>::iterator i=M.begin();
  double num=0;
  while(i!=M.end())
    {
      num=num+1;
      i->value()=num;
      ++i;
    }
  
  std::ofstream out ("sparsity_pattern.1");
  sparsity_pattern.print_gnuplot(out);

  std::ofstream sparsematrix ("sparse_matrix.1");
  A.print(sparsematrix);

  std::ofstream sparsematrix2 ("sparse_matrix.2");
  M.print(sparsematrix2); 

  std::vector<double> x(10,1);
  x=multiply(M,x);
  for(int i=0;i<x.size();i++)
    {
      std::cout<<x[i]<<" ";
    }
 // std::cout<<"\n"<<1e-03<<"\n";

  MinTrace solve;
  
//  solve.min_trace();
  //sparsity_pattern.print();
  return 0;
}
