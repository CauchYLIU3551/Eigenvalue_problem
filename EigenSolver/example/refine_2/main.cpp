/**
 * @file   refine.cpp
 * @author Robert Lie
 * @date   Wed Mar  7 11:43:40 2007
 *
 * @brief  下面的这个程序，是将求解偏微分方程和局部网格加密简单组合在
 *         了一起。我们用的方程是一个含有间断二次系数的椭圆型方程，在
 *         系数间断的位置，解会有一个弱间断。我们瞎算了一个长得象误差
 *         估计的量来做自适应指示子，您可以试试它的效果，;-)
 *
 */

#include <AFEPack/AMGSolver.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/Operator.h>
#include <AFEPack/BilinearOperator.h>
#include <AFEPack/EasyMesh.h>
#include <AFEPack/HGeometry.h>

#define DIM 2

/// 边界条件
double _u_b_(const double * p)
{
  return sin(p[0] + p[1]);
}

/// 右端项
double _f_(const double * p)
{
  return sin(p[0]) + exp(p[1]);
}

/// 二次系数
double A(const AFEPack::Point<DIM>& p)
{
  if (p[0] > p[1]*p[1]) {
    return 4.0;
  } else if (p[0] > sin(p[1])) {
    return 2.0;
  } else {
    return 1.0;
  }
}

/// 二次问题的刚度矩阵，这一段程序我们已经在前面的文章中解释过了。
class Matrix : public StiffMatrix<DIM,double>
{
public:
  Matrix(FEMSpace<double,DIM>& sp) : StiffMatrix<DIM,double>(sp) {};
  virtual void getElementMatrix(const Element<double,DIM> & ele0,
                const Element<double,DIM> & ele1,
                const ActiveElementPairIterator<DIM>::State
                state=ActiveElementPairIterator<DIM>::EQUAL) {
    int n_ele_dof = elementDof0().size();
    double volume = ele0.templateElement().volume();
    const QuadratureInfo<DIM>& quad_info = ele0.findQuadratureInfo(algebricAccuracy());
    std::vector<double> jacobian = ele0.local_to_global_jacobian(quad_info.quadraturePoint());
    int n_quadrature_point = quad_info.n_quadraturePoint();
    std::vector<AFEPack::Point<DIM> > q_point = ele0.local_to_global(quad_info.quadraturePoint());
    std::vector<std::vector<std::vector<double> > > basis_gradient = ele0.basis_function_gradient(q_point);
    std::vector<std::vector<double> > basis_value = ele0.basis_function_value(q_point);
    for (int l = 0;l < n_quadrature_point;l ++) {
      double Jxw = quad_info.weight(l)*jacobian[l]*volume;
      double a_val = A(q_point[l]);
      for (int j = 0;j < n_ele_dof;j ++) {
        for (int k = 0;k < n_ele_dof;k ++) {
          elementMatrix(j,k) += Jxw*a_val*
            innerProduct(basis_gradient[j][l], basis_gradient[k][l]);
        }
      }
    }
  }
};

/// 在当前网格上求解二次方程，然后计算自适应指示子。关于求解方程的部分
/// 和前面的其他程序是一模一样的，但是计算自适应指示子的部分则是完全手
/// 工做的。我想您自己读懂这个部分会比听我讲清楚这个部分要获益更多，我
/// 就特意把注释都抹去了，;-)
void get_indicator(Indicator<DIM>& ind,
                   RegularMesh<DIM>& mesh)
{
  TemplateGeometry<DIM> triangle_template_geometry;
  triangle_template_geometry.readData("triangle.tmp_geo");
  CoordTransform<DIM,DIM> triangle_coord_transform;
  triangle_coord_transform.readData("triangle.crd_trs");
  TemplateDOF<DIM> triangle_template_dof(triangle_template_geometry);
  triangle_template_dof.readData("triangle.1.tmp_dof");
  BasisFunctionAdmin<double,DIM,DIM> triangle_basis_function(triangle_template_dof);
  triangle_basis_function.readData("triangle.1.bas_fun");

  TemplateGeometry<DIM> twin_triangle_template_geometry;
  twin_triangle_template_geometry.readData("twin_triangle.tmp_geo");
  CoordTransform<DIM,DIM> twin_triangle_coord_transform;
  twin_triangle_coord_transform.readData("twin_triangle.crd_trs");
  TemplateDOF<DIM> twin_triangle_template_dof(twin_triangle_template_geometry);
  twin_triangle_template_dof.readData("twin_triangle.1.tmp_dof");
  BasisFunctionAdmin<double,DIM,DIM> twin_triangle_basis_function(twin_triangle_template_dof);
  twin_triangle_basis_function.readData("twin_triangle.1.bas_fun");

  std::vector<TemplateElement<double,DIM,DIM> > template_element(2);
  template_element[0].reinit(triangle_template_geometry,
                             triangle_template_dof,
                             triangle_coord_transform,
                             triangle_basis_function);
  template_element[1].reinit(twin_triangle_template_geometry,
                             twin_triangle_template_dof,
                             twin_triangle_coord_transform,
                             twin_triangle_basis_function);

  FEMSpace<double,DIM> fem_space(mesh, template_element);

  int n_element = mesh.n_geometry(DIM);
  fem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++) {
    if (mesh.geometry(DIM,i).n_vertex() == 3) {
      fem_space.element(i).reinit(fem_space,i,0);
    } else {
      fem_space.element(i).reinit(fem_space,i,1);      
    }
  }
  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  Matrix stiff_matrix(fem_space);
  stiff_matrix.algebricAccuracy() = 2;
  stiff_matrix.build();

  FEMFunction<double,DIM> u_h(fem_space);
  Vector<double> f_h;
  Operator::L2Discretize(&_f_, fem_space, f_h, 3);

  BoundaryFunction<double,DIM> boundary(BoundaryConditionInfo::DIRICHLET,
                                        1, &_u_b_);
  BoundaryConditionAdmin<double,DIM> boundary_admin(fem_space);
  boundary_admin.add(boundary);
  boundary_admin.apply(stiff_matrix, u_h, f_h);

  AMGSolver solver(stiff_matrix);
  solver.solve(u_h, f_h);

  u_h.writeOpenDXData("u_h.dx");

  /// 从这里开始就是计算自适应指示子的部分了，;-)
  u_int n_side = mesh.n_geometry(1);
  std::vector<bool> sflag(n_side, false);
  std::vector<double> sjump(n_side);
  FEMSpace<double,DIM>::ElementIterator
    the_ele = fem_space.beginElement(),
    end_ele = fem_space.endElement();
  for (u_int i = 0;the_ele != end_ele;++ the_ele, ++ i) {
    const GeometryBM& ele_geo = mesh.geometry(DIM,i);
    u_int n_bnd = ele_geo.n_boundary();
    for (u_int j = 0;j < n_bnd;++ j) {
      u_int sid_idx = ele_geo.boundary(j);
      const GeometryBM& sid_geo = mesh.geometry(1, sid_idx);
      const GeometryBM& vtx0 = mesh.geometry(0, sid_geo.vertex(0));
      const AFEPack::Point<DIM>& p0 = mesh.point(vtx0.vertex(0));
      const GeometryBM& vtx1 = mesh.geometry(0, sid_geo.vertex(1));
      const AFEPack::Point<DIM>& p1 = mesh.point(vtx1.vertex(0));
      AFEPack::Point<DIM> p(0.5*(p0[0] + p1[0]), 0.5*(p0[1] + p1[1]));
      std::vector<double> u_h_grad = u_h.gradient(p, *the_ele);

      if (sflag[sid_idx] == false) {
        sjump[sid_idx] = (u_h_grad[0]*(p0[1] - p1[1]) -
                          u_h_grad[1]*(p1[0] - p0[0]));
        sflag[sid_idx] = true;
      } else {
        sjump[sid_idx] -= (u_h_grad[0]*(p0[1] - p1[1]) -
                           u_h_grad[1]*(p1[0] - p0[0]));
        sflag[sid_idx] = false;
      }
    }
  }

  the_ele = fem_space.beginElement();
  for (u_int i = 0;the_ele != end_ele;++ the_ele, ++ i) {
    const GeometryBM& ele_geo = mesh.geometry(DIM,i);
    u_int n_bnd = ele_geo.n_boundary();
    ind[i] = 0.0;
    for (u_int j = 0;j < n_bnd;++ j) {
      u_int sid_idx = ele_geo.boundary(j);
      if (sflag[sid_idx] == true) continue;
      ind[i] += sjump[sid_idx]*sjump[sid_idx];
    }
  }
}


int main(int argc, char * argv[])
{
  /// 声明几何遗传树，并从Easymesh格式的文件中读入数据
  HGeometryTree<DIM> h_tree;
  h_tree.readEasyMesh(argv[1]);

  /// 在背景网格上建立第一个非正则网格，并均匀加密三次
  IrregularMesh<DIM> irregular_mesh(h_tree);
  irregular_mesh.globalRefine(2);

  do {
    /// 对非正则网格做半正则化和正则化
    irregular_mesh.semiregularize();
    irregular_mesh.regularize(false);

    /// 这就是通过正则化得到的网格
    RegularMesh<DIM>& regular_mesh = irregular_mesh.regularMesh();

    /// 将这个网格输出到数据文件中
    regular_mesh.writeOpenDXData("D.dx");
    std::cout << "Press ENTER to continue or CTRL+C to stop ..." << std::flush;
    getchar();

    /// 下面一段计算用于加密的指示子。
    Indicator<DIM> indicator(regular_mesh);
    get_indicator(indicator, regular_mesh);

    /// 下面的几行调用进行自适应的函数，都是套话。
    MeshAdaptor<DIM> mesh_adaptor(irregular_mesh);
    mesh_adaptor.convergenceOrder() = 1.; /// 设置收敛阶为0
    mesh_adaptor.refineStep() = 2; /// 最多允许加密一步
    mesh_adaptor.setIndicator(indicator);
    mesh_adaptor.tolerence() = 1.0e-06; /// 自适应的忍量
    mesh_adaptor.adapt(); /// 完成自适应
  } while (1);    

  return 0;
}

/**
 * end of file
 *
 */
