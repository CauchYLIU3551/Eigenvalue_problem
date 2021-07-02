/**
 * @file   test.cpp
 * @author Chengyu Liu 
 * @date   Tue 08 June 2021
 *
 * @brief  这是个抛物型方程的例子，非常简单，中间就没有解释太多了
 *
 */

#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

#include <base/exceptions.h>
#include <lac/full_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>


#include <AFEPack/EasyMesh.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/BilinearOperator.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/BoundaryCondition.h>

//#include <trace/mintrace.h>
#include <CG/CGSolver.h>
#include <EigenSolver/Miscellaneous.h>
#include <EigenSolver/EigenSolver.h>
#define DIM 3 
#define PI (4.0*atan(1.0)) 

/// 初值和边值的表达式
double _u_(const double * p)
{
  return sin(PI*p[0]) * sin(PI*p[1]) *sin(PI*p[2]);
  //return sin(PI*p[0]) * sin(2*PI*p[1]);
  //return p[0]*exp(p[1]);
}

/// 右端项 In fact, right hand side vector will not be used in Laplacian eigenvalue problems.
double _f_(const double * p)
{
  return 10+5*PI*PI*_u_(p);
  //return p[0]*p[1] + sin(p[1]);
}

double f(const double * p)
{
  return 0;
}

/// Construct the Stiff matrix in the left hand side.
class Stiff_Matrix : public L2InnerProduct<DIM,double>
{
public:
  Stiff_Matrix(FEMSpace<double,DIM>& sp) :
    L2InnerProduct<DIM,double>(sp, sp) {}
  virtual void getElementMatrix(const Element<double,DIM>& e0,
                                const Element<double,DIM>& e1,
                                const ActiveElementPairIterator< DIM >::State s)
  {
    double vol = e0.templateElement().volume();
    u_int acc = algebricAccuracy();
    const QuadratureInfo<DIM>& qi = e0.findQuadratureInfo(acc);
    u_int n_q_pnt = qi.n_quadraturePoint();
    std::vector<double> jac = e0.local_to_global_jacobian(qi.quadraturePoint());
    //AFEPack::Point<DIM> test1;
    //std::vector<AFEPack::Point<DIM>> test2;
    // 
    // Here always get an error! vector<Point<DIM> >;
    // Reason: Because there is a same name class in deal.ii: Point class. While
    // using Point, it might use the dealii.Point defaultly!!!! So you just add 
    // the namespace to control the class will solve the problem;
    std::vector<AFEPack::Point<DIM>> q_pnt = e0.local_to_global(qi.quadraturePoint());
    //std::cout<<"ATTENTION!This is a test flag!!!"<<std::endl;
    std::vector<std::vector<double> > bas_val = e0.basis_function_value(q_pnt);
    std::vector<std::vector<std::vector<double> > > bas_grad = e0.basis_function_gradient(q_pnt);
    u_int n_ele_dof = e0.dof().size();
    for (u_int l = 0;l < n_q_pnt;++ l) {
      double Jxw = vol*qi.weight(l)*jac[l];
      for (u_int i = 0;i < n_ele_dof;++ i) {
        for (u_int j = 0;j < n_ele_dof;++ j) {
          //elementMatrix(i,j) += Jxw*(bas_val[i][l]*bas_val[j][l]/_dt +
	  //                         innerProduct(bas_grad[i][l], bas_grad[j][l]));
	  elementMatrix(i,j) += Jxw*(innerProduct(bas_grad[i][l], bas_grad[j][l]));
        }
      }
    }
  }
};

// Construct the Mass matrix in the right hand side,
class Mass_Matrix : public L2InnerProduct<DIM,double>
{
public:
  Mass_Matrix(FEMSpace<double,DIM>& sp) :
    L2InnerProduct<DIM,double>(sp, sp){}
  virtual void getElementMatrix(const Element<double,DIM>& e0,
                                const Element<double,DIM>& e1,
                                const ActiveElementPairIterator< DIM >::State s)
  {
    double vol = e0.templateElement().volume();
    u_int acc = algebricAccuracy();
    const QuadratureInfo<DIM>& qi = e0.findQuadratureInfo(acc);
    u_int n_q_pnt = qi.n_quadraturePoint();
    std::vector<double> jac = e0.local_to_global_jacobian(qi.quadraturePoint());
    //AFEPack::Point<DIM> test1;
    //std::vector<AFEPack::Point<DIM>> test2;
    // 
    // Here always get an error! vector<Point<DIM> >;
    // Reason: Because there is a same name class in deal.ii: Point class. While
    // using Point, it might use the dealii.Point defaultly!!!! So you just add 
    // the namespace to control the class will solve the problem;
    std::vector<AFEPack::Point<DIM>> q_pnt = e0.local_to_global(qi.quadraturePoint());
    //std::cout<<"ATTENTION!This is a test flag!!!"<<std::endl;
    std::vector<std::vector<double> > bas_val = e0.basis_function_value(q_pnt);
    std::vector<std::vector<std::vector<double> > > bas_grad = e0.basis_function_gradient(q_pnt);
    u_int n_ele_dof = e0.dof().size();
    for (u_int l = 0;l < n_q_pnt;++ l) {
      double Jxw = vol*qi.weight(l)*jac[l];
      for (u_int i = 0;i < n_ele_dof;++ i) {
        for (u_int j = 0;j < n_ele_dof;++ j) {
          /*elementMatrix(i,j) += Jxw*(bas_val[i][l]*bas_val[j][l]/_dt +
	    innerProduct(bas_grad[i][l], bas_grad[j][l]));*/
	  elementMatrix(i,j) += Jxw*(bas_val[i][l]*bas_val[j][l]);
        }
      }
    }
  }
};

void boundary_condition_apply(const FEMSpace<double, DIM>& sp, 
		SparseMatrix<double>& A, 
		BoundaryConditionAdmin<double,DIM> boundary,
		bool preserve_symmetry = true)
{
  u_int n_dof = sp.n_dof();
  const SparsityPattern& spA = A.get_sparsity_pattern();
  const std::size_t * rowstart = spA.get_rowstart_indices();
  const u_int * colnum = spA.get_column_numbers();
  
  for(u_int i=0; i< n_dof; ++ i)
  {
	  int bm = sp.dofInfo(i).boundary_mark;
  }
 
  int numOfbm = 0; 
  for (u_int i = 0; i < n_dof; ++ i){
    int bm = sp.dofInfo(i).boundary_mark;
    if (bm == 0) continue;
    // For this case only 1 for Drichlet boundary condition.
    if (bm != 1) continue; // Attention, For this case specially.
    
    numOfbm += 1;
    for (u_int j = rowstart[i]+1; j < rowstart[i+1]; ++ j){
	    A.global_entry(j) -= A.global_entry(j);
    }
    if (preserve_symmetry) {
      for (u_int j = rowstart[i] + 1;j < rowstart[i + 1];++ j) {
        u_int k = colnum[j];
        const u_int * p = std::find(&colnum[rowstart[k] + 1],
                                    &colnum[rowstart[k + 1]], i);
        if (p != &colnum[rowstart[k+1]]) {
          u_int l = p - &colnum[rowstart[0]];
          A.global_entry(l) -= A.global_entry(l);
        }
      }
    }
  }
}

int main(int argc, char * argv[])
{
  /// 准备网格
  Mesh<DIM> mesh;
  mesh.readData(argv[1]);

  /// 准备参考单元
  /*
  TemplateGeometry<DIM> tmp_geo;
  tmp_geo.readData("triangle.tmp_geo");
  CoordTransform<DIM,DIM> crd_trs;
  crd_trs.readData("triangle.crd_trs");
  TemplateDOF<DIM> tmp_dof(tmp_geo);
  tmp_dof.readData("triangle.2.tmp_dof");
  BasisFunctionAdmin<double,DIM,DIM> bas_fun(tmp_dof);
  bas_fun.readData("triangle.2.bas_fun");
*/
  TemplateGeometry<DIM> tmp_geo;
  tmp_geo.readData("tetrahedron.tmp_geo");
  CoordTransform<DIM,DIM> crd_trs;
  crd_trs.readData("tetrahedron.crd_trs");
  TemplateDOF<DIM> tmp_dof(tmp_geo);
  tmp_dof.readData("tetrahedron.1.tmp_dof");
  BasisFunctionAdmin<double,DIM,DIM> bas_fun(tmp_dof);
  bas_fun.readData("tetrahedron.1.bas_fun");

  std::vector<TemplateElement<double,DIM,DIM> > tmp_ele(1);
  tmp_ele[0].reinit(tmp_geo, tmp_dof, crd_trs, bas_fun);

  /// 定制有限元空间
  FEMSpace<double,DIM> fem_space;
  fem_space.reinit(mesh, tmp_ele);
  u_int n_ele = mesh.n_geometry(DIM);
  fem_space.element().resize(n_ele);
  for (u_int i = 0;i < n_ele;++ i) {
    fem_space.element(i).reinit(fem_space, i, 0);
  }
  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  /// 准备初值
  FEMFunction<double,DIM> u_h(fem_space);
  Operator::L2Interpolate(&_u_, u_h);

  /// 准备边界条件
  Vector<double> rhs;
  Operator::L2Discretize(&f, fem_space, rhs, 4);

  BoundaryFunction<double,DIM> boundary(BoundaryConditionInfo::DIRICHLET,
                                        1,
                                        &_u_);
  BoundaryConditionAdmin<double,DIM> boundary_admin(fem_space);
  boundary_admin.add(boundary);

  // double t;//


  // double dt = 0.01; /// 简单起见，随手取个时间步长算了

  /// 准备线性系统的矩阵
  /*
  Matrix mat(fem_space, dt);
  mat.algebricAccuracy() = 3;
  mat.build();*/


  StiffMatrix<DIM,double> stiff_matrix(fem_space);
  stiff_matrix.algebricAccuracy() = 4;
  stiff_matrix.build();
  boundary_condition_apply(fem_space, stiff_matrix, boundary_admin);

  MassMatrix<DIM,double> mass_matrix(fem_space);
  mass_matrix.algebricAccuracy() = 4;
  mass_matrix.build();
  boundary_condition_apply(fem_space, mass_matrix, boundary_admin);

  
  /*
  /// 准备右端项
  Vector<double> rhs(fem_space.n_dof());
  FEMSpace<double,DIM>::ElementIterator the_ele = fem_space.beginElement();
  FEMSpace<double,DIM>::ElementIterator end_ele = fem_space.endElement();
  for (;the_ele != end_ele;++ the_ele) {
    double vol = the_ele->templateElement().volume();
    const QuadratureInfo<DIM>& qi = the_ele->findQuadratureInfo(3);
    u_int n_q_pnt = qi.n_quadraturePoint();
    std::vector<double> jac = the_ele->local_to_global_jacobian(qi.quadraturePoint());
    std::vector<AFEPack::Point<DIM>> q_pnt = the_ele->local_to_global(qi.quadraturePoint());
    std::vector<std::vector<double> > bas_val = the_ele->basis_function_value(q_pnt);

    /// 当基函数的值已知情况下，可以使用下面的函数来加速
    std::vector<double> u_h_val = u_h.value(bas_val, *the_ele);
    std::vector<std::vector<double> > u_h_grad = u_h.gradient(q_pnt, *the_ele);
    const std::vector<int>& ele_dof = the_ele->dof();
    u_int n_ele_dof = ele_dof.size();
    for (u_int l = 0;l < n_q_pnt;++ l) {
      double Jxw = vol*qi.weight(l)*jac[l];
      double f_val = _f_(q_pnt[l]);
      for (u_int i = 0;i < n_ele_dof;++ i)
	{
	  rhs(ele_dof[i]) += Jxw*bas_val[i][l]*(u_h_val[l]/dt + f_val);
        }
    }
    }*/

  /// 应用边界条件
  //  boundary_admin.apply(mat, u_h, rhs);

  /// 求解线性系统
  /*
  AMGSolver solver;
  solver.lazyReinit(mat);
  solver.solve(u_h, rhs, 1.0e-08, 50);*/
  //TraceSolver soll(stiff_matrix, mass_matrix);
  CGSolver AAAA;
  EigenSolver solver(stiff_matrix, mass_matrix);
  /*
  std::vector<double>rhs(stiff_matrix.m(),0),x(12,0);
  rhs[0]=-0.107343;
  rhs[1]=-0.107920;
  rhs[2]=0.296721;
  rhs[3]=0.136396;
  rhs[4]=-1.038010;
  rhs[5]=-0.456170;
  rhs[6]=0.694243;
  rhs[7]=1.047650;
  rhs[8]=0.488934;
  rhs[9]=0.065605;
  rhs[10]=-0.113350;
  rhs[11]=-0.309275;
  std::vector<double> temp(6,0);
  std::vector<std::vector<double>> tempX(12,temp);
  solver.X=tempX;
  for(int i=0;i<12;i++)
    {
      for(int j=0;j<6;j++)
	{
	  std::cin>>solver.X[i][j];
	}
    }
  solver.get_MX();
  solver.solve(x, rhs, 1.0e-3, 200);
  std::cout<<"The rhs vector before CG is :\n";
  for(int i=0;i<rhs.size();i++)
    {
      std::cout<<rhs[i]<<" ";
    }
  std::cout<<std::endl;
  std::cout<<"The solution of CG\n";
  for(int i=0;i<x.size();i++)
    {
      std::cout<<x[i]<<" ";
    }
  std::cout<<std::endl;
  */
  
  /////////////////////////////
  
  //solver.mintrace(10, 1.0e-3, 200);
  std::vector<double> x(stiff_matrix.m(),1);
  double lambda;
  solver.PowerSolve(x, lambda);
  /*
  std::cout<<"This is the maximal eigenvalues of AX=lambda Mx;\n";
  for (int i=0;i < x.size();i++)
  {
    std::cout<<x[i]<<" ";
  }
  std::cout<<"\n";
*/
  std::vector<double> x2(stiff_matrix.m(),1);
  solver.IPowerSolve(x2, lambda, 100, 1.e-3);

  int eig_num=5;
  std::vector<double> tempxx(eig_num,0), lam_3(eig_num,0);
  std::vector<std::vector<double>> x3(stiff_matrix.m(),tempxx);
  /*
  for(int k=0;k<eig_num;k++)
  {
    x3[k][k]=1;
  }*/

  for (int i=0;i<x3.size();i++)
  {
    for(int j=0;j<x3[0].size();j++)
    {
      x3[i][j] = ((double) rand() / (RAND_MAX));
    }
  }

  solver.BIPowerSolve(x3, lam_3, eig_num, 1000, 1.e-3);
  std::cout<<"The "<<eig_num<<" smallest eigenvalues are:\n";

  for (int k=0;k<eig_num;k++)
  {
    std::cout<<lam_3[k]<<" ";
  }
  std::cout<<"\n";

  std::cout<<"The matrix are:\n";
  show_matrix(x3);


//  std::cout<<"This is the minimal eigenvalues of AX=lambda Mx;\n";
  /*
  for (int i=0;i < x2.size();i++)
  {
    std::cout<<x2[i]<<" ";
  }
  std::cout<<"\n";
  */
  /*for (int i=0;i<solver.lambda.size();i++)
    {
      std::cout<<solver.lambda[i]<<std::endl;
    }
  std::cout<<"\n";
*/
  /*
  std::vector<double> exact_eigen={2.0*PI*PI,PI*PI,PI*PI,0.0};
  std::cout<<" The exact eigens are:\n";
  for (int i=0;i<solver.lambda.size();i++)
    {
      std::cout<<fabs(solver.lambda[i]-exact_eigen[i])<<std::endl;
      // std::cout<<exact_eigen[i]<<std::endl;
    }
  std::cout<<"\n";
  
  std::cout<<"This is the matrix V after the whole process\n";
  for(int i=0;i<solver.X.size();i++)
    {
      for(int j=0;j<solver.X[0].size();j++)
	{
	  std::cout<<solver.X[i][j]<<" ";
	}
      std::cout<<std::endl;
      }*/
  
  
  /*
  std::cout<<"Output the solution matrix X!\n";
  for(int i=0;i<solver.X.size();i++)
    {
      for(int j=0;j<solver.X[0].size();j++)
	{
	  std::cout<<solver.X[i][j]<<" ";
	}
      std::cout<<std::endl;
    }
  std::cout<<std::endl;
  */

  /// 输出数据画图
  /* u_h.writeOpenDXData("u_h.dx");*/
  // std::cout << "Press ENTER to continue or CTRL+C to stop ..." << std::flush;
  // getchar();

  // t += dt; /// 更新时间
    
  //std::cout << "\n\tt = " <<  t << std::endl;

  // Print the stiffness matrix in to .txt and .gnuplot form.
 // std::ofstream sparsematrix2 ("stiff_matrix.1");
 // stiff_matrix.print(sparsematrix2);

  std::filebuf fb;
  fb.open ("stiff_matrix.txt",std::ios::out);
  std::ostream os(&fb);
  stiff_matrix.print_formatted(os, 3, true, 0, "0.0", 1);
  fb.close();


 // Print the mass matrix into .txt and .gnuplot form. 
//  std::ofstream sparsematrix  ("mass_matrix.1");
 // mass_matrix.print(sparsematrix);

  std::filebuf fb2;
  fb2.open ("mass_matrix.txt",std::ios::out);
  std::ostream os2(&fb2);
  mass_matrix.print_formatted(os2, 3, true, 0, "0.0", 1);
  fb2.close();


 /* 
  std::cout<<"Attention! This is print out the columns of the matrix A and M\n";
  std::vector<double> tempx(stiff_matrix.n(),0), tempAx, tempMx;

  for(int i=0;i<solver.A->m();i++)
    {
      tempx[i]=1;
      
      tempAx=multiply(solver.A, tempx);
      for(int j=0;j<tempAx.size();j++)
	{
	  std::cout<<tempAx[j]<<" ";
	}
	std::cout<<std::endl;
      
      tempMx=multiply(solver.M, tempx);
      for(int j=0;j<tempMx.size();j++)
	{
	  std::cout<<tempMx[j]<<" ";
	}
      std::cout<<std::endl;
      tempx[i]=0;
    }
    */
  
  return 0;
}

/**
 * end of file
 *
 */
