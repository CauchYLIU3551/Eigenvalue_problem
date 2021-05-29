/*
 * @file         ELOB.h
 * @author    Chengyu Liu
 * @date       Thur Apr 29 2021
 *
 * @brief       ELOBP4dCG method to get the first p smallest eigenvalues of linear 
 *                  response eigenvalue problem.
 */

#ifndef _elobpcg_h_
#define _elobpcg_h_

#include<math.h>
#include<iostream>
#include<vector>
#include<base/exceptions.h>
#include<lac/full_matrix.h>
#include<lac/sparsity_pattern.h>
#include<lac/sparse_matrix.h>
#include<assert.h>
#include<Mathtools.h>
class elobpcg
{
 public:
  typedef dealii::SparseMatrix<double> Matrix;
  elobpcg();
  elobpcg(const Matrix& k, const Matrix& m);
  elobpcg(const Matrix& k, const Matrix& m, const Matrix& e);

  std::vector<double> proj_K_shift(double shift, std::vector<std::vector<double>> Q, std::vector<double> u, std::vector<double> v); // compute the projection vector of v in to u within K-inner product, with shift and the shift matrix Q;
  void GS_K(double shift, std::vector<std::vector<double>> Q, std::vector<std::vector<double>>& X);

  std::vector<std::vector<double>> initial_Z(std::vector<std::vector<double>> Z);
  const Matrix* K;
  const Matrix* M;
  const Matrix* E;
};

#endif
