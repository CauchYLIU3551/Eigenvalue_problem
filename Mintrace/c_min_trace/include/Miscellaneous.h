#include<cmath>
#include<lac/sparsity_pattern.h>
#include<lac/sparse_matrix.h>
// Multiplication of SparseMatrix A and vector x0;
std::vector<double> multiply(const SparseMatrix<double>* A, std::vector<double> x0)
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
	  SparseMatrix<double>::const_iterator i=A->begin(k);
	  
	  while(i!=A->end(k))
	    {
	      x[k]+=i->value()*x0[i->column()];
	      ++i;
	    }
	}
      return x;
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
