#include <trace/mintrace.h>
#include <CG/CGSolver.h>
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
      i->value()=1;
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
  
  //CGSolver aaa;
  //SparseMatrix<double> A;
  TraceSolver B;
  TraceSolver C(A,M);
  std::vector<double> temp(4,0);
  std::vector<std::vector<double>>a(4,temp),Q(4,temp);
  Q[0][0]=1;
  Q[1][1]=1;
  Q[2][2]=1;
  Q[3][3]=1;
  a[0][0]=4;
  a[0][1]=1;
  a[0][2]=-2;
  a[0][3]=2;
  a[1][0]=1;
  a[1][1]=2;
  a[1][3]=1;
  a[2][0]=-2;
  a[2][2]=3;
  a[2][3]=-2;
  a[3][0]=2;
  a[3][1]=1;
  a[3][2]=-2;
  a[3][3]=-1;
  /*
  for(int i=0;i<a.size();i++)
    {
      for(int j=0;j<a.size();j++)
	{
	  std::cout<<a[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"\n";
  
  B.QRSolver(a);
  for(int i=0;i<a.size();i++)
    {
      for(int j=0;j<a.size();j++)
	{
	  std::cout<<a[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"\n";
  */

  std::vector<double> testx(10,1);
  std::vector<std::vector<double>>V;
  // C.rand_V(1, V);
  /*  std::cout<<"This is C.X\n";
  for(int i=0;i<C.X.size();i++)
    {
      for(int j=0;j<C.X[0].size();j++)
	{
	  std::cout<<C.X[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<std::endl;
  std::cout<<"This is C.MX\n";
  for(int i=0;i<C.X.size();i++)
    {
      for(int j=0;j<C.X[0].size();j++)
	{
	  std::cout<<C.MX[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
    std::cout<<std::endl;*/
  // C.get_MX();
    //C.get_Px(testx);
/*
  for(int j=0;j<testx.size();j++)
    {
      std::cout<<testx[j]<<" ";
    }*/

  
  C.mintrace(10,1.0e-5,20);
  std::cout<<"This is the eigenvalues of AX=lambda Mx;\n";
  for (int i=0;i<C.lambda.size();i++)
    {
      std::cout<<C.lambda[i]<<std::endl;
    }
  std::cout<<"\n";
  


      
  //std::vector<double> tempc(3,0);
  //std::vector<std::vector<double>> tempC(10,tempc);
  /*
  tempC[0][0]=0.27219;
  tempC[0][1]=0.40825;
  tempC[0][2]=-0.87135;
  tempC[1][0]=-0.46837;
  tempC[1][1]=-0.40825;
  tempC[1][2]=-0.33758;
  tempC[2][0]=0.40298;
  tempC[2][1]=-0.40825;
  tempC[2][2]=-0.06539;
  C.X=tempC;
  for(int i=0;i<C.X.size();i++)
    {
      for(int j=0;j<C.X[0].size();j++)
	{
	  std::cout<<C.X[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"\n";

  
  std::cout<<"\n";
  C.get_MX();
  
  for(int i=0;i<C.MX.size();i++)
    {
      for(int j=0;j<C.MX[0].size();j++)
	{
	  std::cout<<C.MX[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"\n";
  
  std::vector<double> RHS(10,0);
  RHS[0]=  -9.9584e-17;
  RHS[1]=    2.0913e-17;
  RHS[2]= -2.1501e-17;
  RHS[3]= -2.6410e-17;
  RHS[4]=3.2734e-17;
  RHS[5]=-9.0672e-04;
  std::vector<double> test(10,1), test2(10,1);
 
  C.get_Ap(test);
  for(int i=0;i<test.size();i++)
    {
      std::cout<<C.Ap[i]<<"\n";
    }
	std::cout<<"\n";
test2=multiply(C.M,test2);
  for(int i=0;i<test.size();i++)
    {
      std::cout<<test2[i]<<"\n";
    }
  */

  
  /*
  B.Householder(a);
  for(int i=0;i<a.size();i++)
    {
      for(int j=0;j<a.size();j++)
	{
	  std::cout<<a[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  //identitymatrix(Q,a.size());
  B.QR(a, Q);
  for(int i=0;i<a.size();i++)
    {
      for(int j=0;j<a.size();j++)
	{
	  std::cout<<a[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<"\n";
    for(int i=0;i<a.size();i++)
    {
      for(int j=0;j<a.size();j++)
	{
	  std::cout<<Q[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  */
 
  std::cout<<"Hello world!\n";
  //  std::ofstream sparsematrix1 ("original_matrix.1");
  // A.print(sparsematrix1);
      // std::ofstream sparsematrix2 ("original_matrix.2");
      // M.print(sparsematrix2);
  return 0;
}
