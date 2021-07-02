#include "CG/CGSolver.h"
#include <EigenSolver/EigenSolver.h>
#include <EigenSolver/Miscellaneous.h>
void fun2(int *a,int **b)
{
  //b=&a;
  **b=*a;
  std::cout<<"This is value of b"<<**b<<std::endl;
}

void fun3(int &w, int **b)
{
  **b=w;
}

int main()
{
  int a=2312321,c=2;
  int *b=&c;
  //b=&a;
  // fun3(a,&b);
  std::cout<<*b<<std::endl;
  CGSolver AAA;
  dealii::SparseMatrix<double> A, A2, M, M2;
  // std::ofstream sparsematrix1 ("original_matrix.1");
  //A.print(sparsematrix1);
  std::vector<unsigned int> row_length(10,10);
  unsigned int t=3;
  dealii::SparsityPattern sparsity_pattern(10,10,{2,3,3,3,3,3,3,3,3,2});
  dealii::SparsityPattern sp2(3,3,{2,3,2});
  dealii::SparsityPattern sp3(3,3,{1,1,1}), sp4(10,10,{1,1,1,1,1,1,1,1,1,1});
  sp2.add(0,1);
  sp2.add(1,0);
  sp2.add(1,2);
  sp2.add(2,1);
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
  A2.reinit(sp3);
  M2.reinit(sp4);
    // In this way, I can construct a SparseMatrix to be the test data for the algorithm.
  for (int k=0;k<A.m();k++)
    {
      dealii::SparseMatrix<double>::iterator i=A.begin(k);
      i->value()=2;
      while(i!=A.end(k))
	{
	  i->value()+=1;
	  ++i;
	}
    }

      // In this way, I can construct a SparseMatrix to be the test data for the algorithm.

 // To construct the matrix M as diag(2,2,2) with subdiagonal entries -1. 
 /*
  for (int k=0;k<M.m();k++)
    {
      dealii::SparseMatrix<double>::iterator i=M.begin(k);
      i->value()=3;
      while(i!=M.end(k))
	{
	  i->value()-=1;
	  ++i;
	}
    }
    */
  for (int k=0;k<M.m();k++)
    {
      dealii::SparseMatrix<double>::iterator i=M.begin(k);
      i->value()=20.0*(k+1);
      while(i!=M.end(k))
        {
          i->value()-=1;
	  ++i;
        }
    }

  for (int k=0;k<M2.m();k++)
    {
      dealii::SparseMatrix<double>::iterator i=M2.begin(k);
      i->value()=k+1;
      
    }


 // Construct the right hand side matrix in example as A2. 
    for (int k=0;k<A2.m();k++)
    {
      dealii::SparseMatrix<double>::iterator i=A2.begin(k);
      i->value()=(4*k)+1;
      //i->value()=1;
    }



  //std::cout<<AAA.tolerence()<<std::endl;
  //CGSolver solve(A), sol2(M);
  EigenSolver sol3(M,A2),sol4(A,M2);
  std::ofstream sparsematrix2 ("sparse_matrix.2");
  A.print(sparsematrix2);
  //std::cout<<solve.A->n()<<std::endl;
  dealii::Vector<double> x0(10),B(10), x1(3),B2(3);
  for(int i=0;i<x1.size();i++)
  {
    std::cout<<x1[i]<<" ";
  }
  std::cout<<"\n";
  B[0]=4;
  for(int i=1;i<9;i++)
    {
      B[i]=5;
    }
  B[9]=4;
  //B2[0]=1;
  //B2[2]=1.8;
  B2[0]=-2;
  B2[1]=4;
  B2[2]=-2;
  std::cout<<"test1\n";
  //solve.solve(x0,B,1.0e-12,100);
  //sol2.solve(x1,B2,1.0e-3,20);
  /*
  for (int k=0;k<10;k++)
    {
      std::cout<<x0[k]<<" ";
    }
  std::cout<<"\n";
  
  for (int k=0;k<3;k++)
    {
        std::cout<<x1[k]<<" ";
    }
  std::cout<<"\n";
  */
  std::vector<double> xx(3,1);
  double lambda = 0;
  std::cout<<"11\n";
  sol3.PowerSolve(xx, lambda);
  std::cout<<"22\n";
  for (int k=0;k<xx.size();k++)
  {
    std::cout<<xx[k]<<" ";
  }
  std::cout<<"\n";
  std::vector<double> x2(3,1);
  sol3.IPowerSolve(x2,lambda);
    for (int k=0;k<x2.size();k++)
  {
    std::cout<<x2[k]<<" ";
  }
  std::cout<<"\n";
  int eig_num=3;
  std::vector<double> tempxx(eig_num, 0), lam_3(eig_num,0);
  std::vector<std::vector<double>> x3(3, tempxx);
  x3[0][0]=1;
  x3[1][1]=1;
  x3[2][2]=1;
  for (int i=0;i<x3.size();i++)
  {
    for(int j=0;j<x3[0].size();j++)
    {
      x3[i][j] = ((double) rand() / (RAND_MAX));
    }
  }
  std::cout<<"The initial matrix x3 is :: \n";
  show_matrix(x3);
  sol3.BIPowerSolve(x3,lam_3,eig_num,100);
    for (int k=0;k<eig_num;k++)
  {
    std::cout<<lam_3[k]<<" ";
  }
  std::cout<<"\n";
  std::cout<<"The eigenvectors are:\n";
  show_matrix(x3);




  eig_num=3;
  std::vector<double> tempx4(eig_num, 0), lam_4(eig_num,0);
  std::vector<std::vector<double>> x4(A.m(), tempxx);
  //x3[0][0]=1;
  //x3[1][1]=1;
  //x3[2][2]=1;
  for (int i=0;i<x4.size();i++)
  {
    for(int j=0;j<x4[0].size();j++)
    {
      x4[i][j] = ((double) rand() / (RAND_MAX));
    }
  }
  std::cout<<"The initial matrix x4 is :: \n";
  show_matrix(x4);
  sol4.LOBPCGSolve(x4,lam_4,eig_num,100);
    for (int k=0;k<eig_num;k++)
  {
    std::cout<<lam_4[k]<<" ";
  }
  std::cout<<"\n";



 // std::cout<<"The eigenvectors are:\n";
 // show_matrix(x3);
 

  // std::cout<<(*P).n()<<std::endl;
  std::ofstream sparsematrix ("sparse_matrix.2");
  M.print(sparsematrix);
  std::ofstream sparsematrixA ("sparse_matrixA");
  A2.print(sparsematrixA);

  std::ofstream sparsematrix4 ("matrixA");
  A.print(sparsematrix4);


  std::filebuf fb;
  fb.open ("stiff_matrix.txt",std::ios::out);
  std::ostream os(&fb);
  A.print_formatted(os, 3, true, 0, "0.0", 1);
  fb.close();

  std::filebuf fb2;
  fb2.open ("mass_matrix.txt",std::ios::out);
  std::ostream os2(&fb2);
  M2.print_formatted(os2, 3, true, 0, "0.0", 1);
  fb2.close();




  return 0;
}
