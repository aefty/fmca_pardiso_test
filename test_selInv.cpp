// This file is part of FMCA, the Fast Multiresolution Covariance Analysis
// package.
//
// Copyright (c) 2022, Michael Multerer
//
// All rights reserved.
//
// This source code is subject to the GNU Affero General Public License v3.0
// license and without any warranty, see <https://github.com/muchip/FMCA>
// for further information.
//
// #define EIGEN_DONT_PARALLELIZE
#include <FMCA/src/util/Tictoc.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <FMCA/CovarianceKernel>
#include <FMCA/Samplets>
#include <iostream>
#include <iomanip>

#include "pardiso_interface.h"

//#define NPTS 100000
#define DIM 2

using Interpolator = FMCA::TotalDegreeInterpolator;
using SampletInterpolator = FMCA::MonomialInterpolator;
using Moments = FMCA::NystromMoments<Interpolator>;
using SampletMoments = FMCA::NystromSampletMoments<SampletInterpolator>;
using MatrixEvaluator = FMCA::NystromEvaluator<Moments, FMCA::CovarianceKernel>;
using H2SampletTree = FMCA::H2SampletTree<FMCA::ClusterTree>;
using SparseMatrix = Eigen::SparseMatrix<FMCA::Scalar, Eigen::RowMajor>;



int main(int argc, char *argv[]) {

  int NPTS    = std::atoi(argv[1]);
  int k       = std::atoi(argv[2]);
  double reg  = std::atof(argv[3]);
  double r    = std::atof(argv[4]);
  int verbose = std::atoi(argv[5]);

  FMCA::Tictoc T;
  const FMCA::Index dtilde = 4;
  const FMCA::Index mpole_deg = 2 * (dtilde - 1);
  const FMCA::Scalar eta = 0.5;
  const FMCA::CovarianceKernel function("EXPONENTIAL", r);
  const FMCA::Matrix P = 0.5 * (FMCA::Matrix::Random(DIM, NPTS).array() + 1);
  const FMCA::Scalar threshold = 1e-12;
  const Moments mom(P, mpole_deg);
  const MatrixEvaluator mat_eval(mom, function);

  std::cout << "npts:                         " << NPTS << std::endl;
  std::cout << "dim:                          " << DIM << std::endl;
  std::cout << "dtilde:                       " << dtilde << std::endl;
  std::cout << "mpole_deg:                    " << mpole_deg << std::endl;
  std::cout << "eta:                          " << eta << std::endl;
  const SampletMoments samp_mom(P, dtilde - 1);
  H2SampletTree hst(mom, samp_mom, 0, P);
  T.tic();
  FMCA::internal::SampletMatrixCompressor<H2SampletTree> Scomp;
  Scomp.init(hst, eta, threshold);
  T.toc("planner:                     ");
  T.tic();
  Scomp.compress(mat_eval);
  T.toc("compressor:                  ");
  T.tic();
  const auto &trips = Scomp.triplets();
  T.toc("triplets:                    ");
  std::cout << "anz:                          "
            << std::round(trips.size() / FMCA::Scalar(NPTS)) << std::endl;
  
  SparseMatrix S(NPTS, NPTS);
  S.setFromTriplets(trips.begin(), trips.end());
  FMCA::Vector x(NPTS), y1(NPTS), y2(NPTS);
  FMCA::Scalar err = 0;
  FMCA::Scalar nrm = 0;
  for (auto i = 0; i < 10; ++i) {
    FMCA::Index index = rand() % P.cols();
    x.setZero();
    x(index) = 1;
    FMCA::Vector col = function.eval(P, P.col(hst.indices()[index]));
    y1 = col(Eigen::Map<const FMCA::iVector>(hst.indices(), hst.block_size()));
    x = hst.sampletTransform(x);
    y2.setZero();
    for (const auto &i : trips) {
      y2(i.row()) += i.value() * x(i.col());
      if (i.row() != i.col()) y2(i.col()) += i.value() * x(i.row());
    }
    FMCA::Vector y3 = S.selfadjointView<Eigen::Upper>() * x;
    assert((y3 - y2).norm() / y2.norm() && "this needs to be small");
    y2 = hst.inverseSampletTransform(y2);
    err += (y1 - y2).squaredNorm();
    nrm += y1.squaredNorm();
  }
  err = sqrt(err / nrm);
  std::cout << "compression error:            " << err << std::endl
            << std::flush;
  std::cout << std::string(72, '-') << std::endl;
  SparseMatrix I(NPTS, NPTS);
  S.setFromTriplets(trips.begin(), trips.end());
  // add regularization (Pardiso might be sensitive to this)
  I.setIdentity();
  S += I*reg;
  S.makeCompressed();



  //S.setIdentity();

  Eigen::SparseMatrix<FMCA::Scalar, Eigen::RowMajor> invS = S;
  // better save than sorry
  invS.makeCompressed();
  std::cout << "Eigen::Sparse written         " << std::endl;
  std::cout << "invS norm: " << invS.norm() << std::endl;
  SparseMatrix::StorageIndex n = invS.rows();
  SparseMatrix::StorageIndex *ia = invS.outerIndexPtr();
  SparseMatrix::StorageIndex *ja = invS.innerIndexPtr();
  SparseMatrix::Scalar *a = invS.valuePtr();
  std::cout << "entering pardiso block\n" << std::flush;
  T.tic();
  std::printf("ia=%p ja=%p a=%p n=%i nnz=%i\n", ia, ja, a, n, ia[n]);
  std::cout << std::flush;
  pardiso_interface(ia, ja, a, n);
  std::printf("ia=%p ja=%p a=%p n=%i nnz=%i\n", ia, ja, a, n, ia[n]);
  std::cout << std::string(72, '-') << std::endl;
  T.toc("Wall time pardiso:           ");

  FMCA::Matrix X(NPTS, k);
  X.setRandom();

  FMCA::Matrix Y = invS.selfadjointView<Eigen::Upper>() * (S.selfadjointView<Eigen::Upper>() * X).eval();
  std::cout << "error: " << (X - Y).norm() / X.norm() << std::endl;
  
  Eigen::SparseMatrix<FMCA::Scalar> invS_sym = Eigen::SparseMatrix<FMCA::Scalar>(invS.selfadjointView<Eigen::Upper>());
  Eigen::SparseMatrix<FMCA::Scalar> S_sym    = Eigen::SparseMatrix<FMCA::Scalar>(S.selfadjointView<Eigen::Upper>());
  Eigen::SparseMatrix<FMCA::Scalar> temp     = invS_sym * S_sym;
  double norm = temp.norm();
  std::cout << "error inv : " << std::abs(norm*norm-  temp.rows())/temp.rows()  << std::endl;


  if(verbose>0){
  // Create a view of the top-left 10x10 submatrix and convert it to dense
  int b = 10;
  Eigen::MatrixXd invS_exact_dense = S_sym.toDense().inverse();
  Eigen::MatrixXd invS_exact_dense_temp = invS_exact_dense.block(0, 0, b, b);
  std::cout <<"---- invS Exact*100 ---- \n"  << std::fixed << std::setprecision(1)<< 100*invS_exact_dense_temp << std::endl;

  Eigen::MatrixXd invS_temp = invS_sym.block(0, 0, b, b).toDense();
  std::cout <<"---- invS*100 ----\n"  << std::fixed << std::setprecision(1) << 100*invS_temp << std::endl;

  Eigen::MatrixXd S_temp = S_sym.block(0, 0, b, b).toDense();
  std::cout <<"---- S ---- \n"  << std::fixed << std::setprecision(1)<< S_temp << std::endl;

  }

  return 0;
}
