#include <iostream>
#include <fstream>
#include <iterator>
#include <mintrace.h>
// Following header files are used in the step-3 in the examples of deal.II;
// Here just using these to get a large scale matrix as the test data.
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
// Finally, this is for output to a file and to the console:
#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>

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

// This function computes the inner-product corresponding with M of <u,v>
void M_inner(SparseMatrix<double> M, std::vector<double>u, std::vector<double>v){}


class Step3
{
public:
  Step3 ();

  void run ();
  void make_grid ();
  void setup_system ();
  void assemble_system ();

  //void multiply(std::vector<double>& x);
  std::vector<double> multiply(std::vector<double> x0);

  
  //void solve ();
  //
  void output_results () const;

  Triangulation<2>     triangulation;
  FE_Q<2>              fe;
  DoFHandler<2>        dof_handler;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       system_rhs;
  
  //std::vector<double> x;

};



//
//
//void Step3::multiply(std::vector<double>& x)
std::vector<double> multiply(SparseMatrix<double>A, std::vector<double> x0)
{
  std::vector<double> x(x0.size(),0);
  
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



Step3::Step3 ()
  :
  fe (1),
  dof_handler (triangulation)
{}

void Step3::make_grid ()
{
  // First create the grid and refine all cells five times. Since the initial
  // grid (which is the square [-1,1]x[-1,1]) consists of only one cell, the
  // final grid has 32 times 32 cells, for a total of 1024.
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (5);
  // Unsure that 1024 is the correct number?  Let's see: n_active_cells
  // returns the number of active cells:
  std::cout << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl;
  
  std::cout << "Total number of cells: "
            << triangulation.n_cells()
            << std::endl;
  // Note the distinction between n_active_cells() and n_cells().
}

void Step3::setup_system ()
{
  dof_handler.distribute_dofs (fe);
  std::cout << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << std::endl;
  
  CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, c_sparsity);
  sparsity_pattern.copy_from(c_sparsity);

  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}

void Step3::assemble_system ()
{

  QGauss<2>  quadrature_formula(2);

  FEValues<2> fe_values (fe, quadrature_formula,
                         update_values | update_gradients | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  DoFHandler<2>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

      cell_matrix = 0;
      cell_rhs = 0;

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            cell_matrix(i,j) += (fe_values.shape_grad (i, q_point) *
                                 fe_values.shape_grad (j, q_point) *
                                 fe_values.JxW (q_point));

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          cell_rhs(i) += (fe_values.shape_value (i, q_point) *
                          1 *
                          fe_values.JxW (q_point));

      cell->get_dof_indices (local_dof_indices);
      ///////////////////////
      // Following commands compute the interior of the matrix and store them into the matrix.
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             cell_matrix(i,j));

      // And again, we do the same thing for the right hand side vector.
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }

  std::map<types::global_dof_index,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<2>(),
                                            boundary_values);

  MatrixTools::apply_boundary_values (boundary_values,
                                      system_matrix,
                                      solution,
                                      system_rhs);
}

//void Step3::multiply(std::vector<double>& x)
std::vector<double> Step3::multiply(std::vector<double> x0)
{
  std::vector<double> x(x0.size(),0);
  
  if (system_matrix.n()!=x0.size())
    {
      std::cout<<"Function multiply() ERROR: The sizes of matrix and vector are not same! Please check it!\n";
    }
  else
    {
      // std::cout<<"TEST POINT 2 \n";
      for(int k=0;k<system_matrix.m();k++)
	{
	  //std::cout<<"this is the "<<k<<"th iterations \n";
	  SparseMatrix<double>::const_iterator i=system_matrix.begin(k);
	  
	  while(i!=system_matrix.end(k))
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



void Step3::run ()
{
  make_grid ();
  setup_system();
  assemble_system ();
  std::cout<<system_matrix.m()<<"\n";
  std::cout<<system_matrix.n()<<"\n";

  //std::cout<<system_matrix.begin()->value()<<"\n";
  //std::cout<<system_matrix.begin()->column()<<"\n";

  std::vector<double> x(system_matrix.n(),1);
  //
  // Due to the system_matrix is a member of the class step3, but x is not a member of the class,
  // So we can not directly define the multiply function out of the class and using the
  // system_matrix directly.
  // we need to get a variable A equals to the system_matrix but A is not a member of the class;
  
  std::cout << "CheckPoint 1 \n";

  // multiply(x);
  x=multiply(x);
  //std::cout<<"This is the value of the vector x"<<x[0]<<std::endl;
  std::ofstream out ("sparse_matrix");
  system_matrix.print(out);

  std::ofstream output_file ("solution");
  for (const auto &e : x) output_file << e << "\n";
  
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

// The function using Conjugate Gradient method to get the solution of the problem Ax=b;
// It is modified CG function that is designed for solve PAP delta = PAX which needs A, M, X, b;
// ATTENTION: It can be revised as a class for A and M, It will be more convenient.
std::vector<double> CG(SparseMatrix<double>A, SparseMatrix<double> M, std::vector<std::vector<double>> X, std::vector<double> b, std::vector<double> &x)
{
  return x;
}

// Computing the matrix P during the minimal trace process;
std::vector<std::vector<double>> matrixP(SparseMatrix<double> M, std::vector<std::vector<double>>X )
{
  std::vector<std::vector<double>> P;
  return P;
}

//
// This function computes the matrix M-orthogonal modified Gram-Schmidt procedure;
std::vector<double> M_orth(SparseMatrix<double> M, std::vector<double> x)
{
  return x;
}

// This function computes a random matrix V as the initial matrix of the iteration;
// n implies the column number of V;
void rand_V(SparseMatrix<double>M, std::vector<std::vector<double>> &V, int n)
{
  std::vector<std::vector<double>> A(V);
  V=A;
}

// Compute matrix A - B;
std::vector<std::vector<double>> minus(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B)
{
  return A;
}

// This function computes the smallest p eigenvalues of matrix A corresponding to matrix M; 
// Notice: The A and M are square matrices here.
std::vector<double> min_trace(SparseMatrix<double> A, SparseMatrix<double>M, int p, double tol=0.001, int max_iter=1)
{
  int n=A.m();
  std::vector<std::vector<double>>V(p);
  std::vector<std::vector<double>>W(n);
  std::vector<std::vector<double>> MXTheta(p), U, X, Rk, delta;
  //SparseMatrix<double> P;
  
  // Construct the initial matrix V1 which is a n x p matrix and it is orthogonal by M, i.e.
  // V1'*M*V1=I_p;
  // I store the V(i,j)=V[j][i];
  rand_V(M,V,p);
  for (int k=0;k<max_iter;k++)
    {
      // Computing the matrix W and H;
      // Here I use the property of the symmetric A: H=V'*A*V=V'*W=W'*V=H'=V'*A'*V=H;

      // store the A *V[i] into W[i];
      for(int i=0;i<p;i++)
	{
	  std::vector<double> tempV=V[i];
	  // W[i]=multiply(A,V[i]);
	  W[i]=multiply(A,tempV);
	}

      // copy the W into H;
      std::vector<std::vector<double>> H(W);
      // Computing H'*V; in fact H[i][j]=H(j,i); So above the H=H';
      // Get the actual Hk;
      multiply(H,transpose(V));

      // Compute the spectral decomposition of Hk;
      identitymatrix(U,p);
      QRSolver(H,U); // H = theta_k, U=U_k;

      // transpose V into V[i][j]=V(i,j);
      X=transpose(V);
      // Compute Ritz vectors X_k=V_k*U_k;
      multiply(X,U);
      
      // Compute MXTheta = M * X_k * Theta_k
      // using XTheta to represent the matrix otherwise X will be changed by multiply;
      std::vector<std::vector<double>> XTheta(X);    
      // Compute XTheta = X_k * Theta_k = XTheta*H;
      multiply(XTheta,H);
      for(int i=0;i<p;i++)
	{
	  MXTheta[i]=multiply(M, XTheta[i]);
	}
      // Compute residual Rk in fact Rk[i][j] = Rk(j,i) ;
      // This step might be wrong;
      //Rk=minus(multiply(transpose(U),W), MXTheta);
      /*
      if(infi_Norm(Rk)<tol)
	{
	  break;
	}
      */

      // Compute matrix P;
      //  P=matrixP(M, X);
      
      // Solve the SPD eigenvalue problem by modified PCG or CG;
      // PAP delta_k = PA X_k;
      // divide it into p eigenvalue problem:
      // PAP d_i = PA x_i; d_i = delta_k*ei, x_i = X_k*e_i;
      
      for(int i=0;i<p;i++)
	{
	  std::vector<double> x(n,0);
	  SparseMatrix<double> PAP;
	  std::vector<double> PAXi;
	  delta[i]=CG(A, M, X, PAXi, x);
	}


      // update the V = V_k+1;
      std::vector<std::vector<double>> XT;
      XT=transpose(X);
      for(int i=0;i<p;i++)
	{
	  XT=minus(XT, delta);
	  V[i]=M_orth(M,XT[i]);
	}
      //V=transpose(V);
    }

  // Compute the diagonal Matrix V'*A*V whose diagonal entries are eigenvalues of the matrix A
  // corresponding to matrix M;

  std::vector<double> eigenvalue(p);
  for(int i=0;i<p;i++)
    {
      std::vector<double> tempV, tempV2;
      tempV=multiply(A,V[i]);
      tempV2=V[i];
      eigenvalue[i]=multiply(tempV,tempV2);
    }

  return eigenvalue;

}



int main()
{
  Step3 laplace_problem;
  /////laplace_problem.run ();

  
  //std::stringstream result;
  std::vector<double> b(4,0);
  std::vector<std::vector<double>> A(4, b), P(4,b);
  P[0][0]=1;
  P[1][1]=1;
  P[2][2]=1;
  P[3][3]=1;
  std::vector<std::vector<double>> Q(P);

  A[0][0]=4;
  A[0][1]=1;
  A[0][2]=-2;
  A[0][3]=2;
  A[1][0]=1;
  A[1][1]=2;
  A[1][2]=0;
  A[1][3]=1;
  A[2][0]=-2;
  A[2][1]=0;
  A[2][2]=3;
  A[2][3]=-2;
  A[3][0]=2;
  A[3][1]=1;
  A[3][2]=-2;
  A[3][3]=-1;

  
  //std::vector<double> c={1,2,3,4};
  //transpose(A);

  std::cout<<"The matrix A is :\n";
  std::cout<<"\n";
 
      
  for(int i=0;i<A.size();i++)
    {
      for(int j=0; j<A.size();j++)
	{
	  std::cout<<A[i][j]<<" ";
	}
      std::cout<<"\n";
    }
  std::cout<<"\n";
  
  // multiply(A,A);
  
  // c=multiply(A,c);

  /*
  Householder(A,P);

  std::cout<<"The matrix A is :\n";
  std::cout<<"\n";
      
  for(int i=0;i<A.size();i++)
    {
      for(int j=0; j<A.size();j++)
	{
	  std::cout<<A[i][j]<<" ";
	}
      std::cout<<"\n";
    }
  std::cout<<"\n";

  QR(A,Q);

  */

  ////
  //  ATTENTION:!!!!! the original matrix A is symmetric, but after A= QR to get A_=RQ, A_ might be  // nonsymmetric!!!!!! the function is useless!!!!!!! Please edit the function.
  QRSolver(A,Q);

  std::cout<<"The matrix A is :\n";
  std::cout<<"\n";
      
  for(int i=0;i<A.size();i++)
    {
      for(int j=0; j<A.size();j++)
	{
	  std::cout<<A[i][j]<<" ";
	}
      std::cout<<"\n";
    }
  std::cout<<"\n";

  std::cout<<"The matrix Q is :::::::\n";
  std::cout<<"\n";
      
  for(int i=0;i<Q.size();i++)
    {
      for(int j=0; j<Q.size();j++)
	{
	  std::cout<<Q[i][j]<<" ";
	}
      std::cout<<"\n";
    }
  std::cout<<"\n";

  std::cout<<"check if the matrix Q is orthogonal\n";
  multiply(Q,transpose(Q));

  std::cout<<"The matrix Q*Q' is :::::::\n";
  std::cout<<"\n";
      
  for(int i=0;i<Q.size();i++)
    {
      for(int j=0; j<Q.size();j++)
	{
	  std::cout<<Q[i][j]<<" ";
	}
      std::cout<<"\n";
    }
  std::cout<<"\n";
  
  std::cout<<"hello world!"<<std::endl;
  //SparseMatrix<double> A;
  return 0;
}
