#include<iostream>
#include<cmath>
#include<vector>

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

double infi_norm(std::vector<double> a)
{
  double temp=0;
  for (int i=0;i<a.size();i++)
    {
      if(temp<fabs(a[i]))
	{
	  temp=fabs(a[i]);
	}
    }
  return temp;
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
  //p=res;
  for(int i=0;i<n;i++)
    {
      p[i]=-res[i];
    }
  //Ap=get_Ap(A,p);
   if (infi_norm(p)<tol)
    {
      return;
    }
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
  std::cout<<"the first delta is : "<<delta<<std::endl;
  
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
  while(infi_norm(res)>tol&& iter < max_iter)
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
      std::cout<<"This is "<<iter<<" iterations!\n"
	       <<"x is \n";
      for(int i=0;i<x.size();i++)
	{
	  std::cout<<x[i]<<" ";
	}
      std::cout<<"\n";
    }
}

void fun1(std::vector<double >& x)
{
  x[0]+=1;
}
void fun2(std::vector<double> x)
{
  fun1(x);
  std::cout<<"fun2 !!!!\n";
  for(int i=0;i<x.size();i++)
    {
      std::cout<<x[i]<<" ";
    }
}

int main()
{
  double tol = 1.0e-5;
  // const int n = 3;
  std::vector<double> x(3, 0), rhs(3, 0);
  std::vector<std::vector<double>> A(3,x);
  A[0][0]=2;
  A[0][1]=-1;
  A[1][0]=-1;
  A[1][1]=2;
  A[1][2]=-1;
  A[2][1]=-1;
  A[2][2]=2;
  rhs[0]=1;
  rhs[2]=1.8;
  CG(A, x, rhs, tol);
  std::vector<double> b;
  b=x;
  for(int i=0;i<b.size();i++)
    {
      std::cout<<b[i]<<std::endl;
    }

  std::cout<<" This is a special case to compute CG\n";
  std::vector<double> rhs2(1,0),x2(1,0);
  std::vector<std::vector<double>> AA(1,rhs2);
  CG(AA, x2, rhs2, tol);

  for(int i=0;i<x2.size();i++)
    {
      std::cout<<x2[i]<<std::endl;
      }

  std::vector<double> t(2,0), t2(2,0);
  fun1(t);
  for(int i=0;i<t.size();i++)
    {
      std::cout<<t[i]<<" ";
    }
  fun2(t2);
  
  return 0;
}
