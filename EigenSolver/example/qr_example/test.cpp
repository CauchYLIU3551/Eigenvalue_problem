#include<iostream>
#include<EigenSolver/EigenSolver.h>
#include<EigenSolver/Miscellaneous.h>
int main()
{
  std::vector<double> tempx(5,0);
  std::vector<std::vector<double>> A(5,tempx), X(5,tempx);
  for (int i=0; i< 5;i++)
    {
      X[i][i]=1;
    }
  A[0][0]=75.27960;
  A[0][4]=-56.95960;
  A[4][0]=A[0][4];
  A[1][1]=139.20000;
  A[2][2]=249.67700;
  A[3][3]=153.60000;
  A[4][4]=129.34800;
  std::cout<<"The matrix A is :\n";
  show_matrix(A);
  std::cout<<"The matrix X is :\n";
  show_matrix(X);
  QRSolver(A,X);
  std::cout<<"The matrix A and X after QRSolver are : \n";
  show_matrix(A);
  std::cout<<"The matrix X is :\n";
  show_matrix(X);
  std::cout<<"Hello World!\n";
  return 0;
}

