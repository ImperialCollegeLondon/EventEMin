#ifndef EVENT_EMIN_DATA_STATS_H
#define EVENT_EMIN_DATA_STATS_H

#include <Eigen/SVD>
#include <cassert>

#include "EventEMin/types_def.h"

namespace EventEMin
{
template <typename Derived1, typename Derived2>
void
computeWhitening(const MatrixBase<Derived1>& cov, MatrixBase<Derived2>& w)
{
  const Eigen::JacobiSVD<Derived1> svd(cov, Eigen::ComputeFullU);
  w.noalias() = svd.matrixU() *
                (svd.singularValues().array() + 1.0e-8)
                    .sqrt()
                    .inverse()
                    .matrix()
                    .asDiagonal() *
                svd.matrixU().transpose();
}
template <typename Derived1, typename Derived2, typename Derived3>
void
computeWhitening(const MatrixBase<Derived1>& cov, MatrixBase<Derived2>& w,
                 MatrixBase<Derived3>& sValues)
{
  const Eigen::JacobiSVD<Derived1> svd(cov, Eigen::ComputeFullU);
  sValues = (svd.singularValues().array() + 1.0e-8).sqrt();
  w.noalias() = svd.matrixU() *
                sValues.array().inverse().matrix().asDiagonal() *
                svd.matrixU().transpose();
}
template <typename Derived1, typename Derived2, typename Derived3>
void
colourPoints(const MatrixBase<Derived1>& cWhitened,
             const MatrixBase<Derived2>& w, MatrixBase<Derived3>& cCentred)
{
  cCentred.noalias() = w.inverse() * cWhitened;
}
template <typename Derived1, typename Derived2, typename Derived3>
void
whitenPoints(const MatrixBase<Derived1>& cCentred,
             const MatrixBase<Derived2>& w, MatrixBase<Derived3>& cWhitened)
{
  cWhitened.noalias() = w * cCentred;
}

template <typename T>
class DataStats
{
 protected:
  Vector<T> min_, max_;
  Vector<T> mean_;
  Matrix<T> centred_, cov_;

 public:
  DataStats(void) = default;

  const Vector<T>&
  min(void) const
  {
    return min_;
  }
  const Vector<T>&
  max(void) const
  {
    return max_;
  }
  const Vector<T>&
  mean(void) const
  {
    return mean_;
  }
  const Matrix<T>&
  centred(void) const
  {
    return centred_;
  }
  const Matrix<T>&
  cov(void) const
  {
    return cov_;
  }

  void
  computeAll(const Ref<const Matrix<T> >& data)
  {
    computeLimits(data);
    computeMoments(data);
  }
  void
  computeAll(const Ref<const Matrix<T> >& data,
             const Ref<const Vector<T> >& weights)
  {
    assert(data.cols() == weights.size());
    computeLimits(data);
    computeMoments(data, weights, weights.sum() + T(1.0e-32));
  }

  void
  computeCentred(const Ref<const Matrix<T> >& data)
  {
    centred_ = data.colwise() - mean();
  }
  void
  computeCov(void)
  {
    assert(centred().cols() > 1);
    cov_.noalias() = centred() * centred().transpose() / (centred().cols() - 1);
  }
  void
  computeCov(const Ref<const Vector<T> >& weights, const T& weightsSum)
  {
    cov_.noalias() =
        (centred().array().rowwise() * weights.transpose().array()).matrix() *
        centred().transpose() / weightsSum;
  }
  void
  computeMean(const Ref<const Matrix<T> >& data)
  {
    mean_ = data.rowwise().mean();
  }
  void
  computeMean(const Ref<const Matrix<T> >& data,
              const Ref<const Vector<T> >& weights, const T& weightsSum)
  {
    mean_.noalias() = (data.array().rowwise() * weights.transpose().array())
                          .matrix()
                          .rowwise()
                          .sum() /
                      weightsSum;
  }
  void
  computeMoments(const Ref<const Matrix<T> >& data)
  {
    computeMean(data);
    computeCentred(data);
    computeCov();
  }
  void
  computeMoments(const Ref<const Matrix<T> >& data,
                 const Ref<const Vector<T> >& weights, const T& weightsSum)
  {
    computeMean(data, weights, weightsSum);
    computeCentred(data);
    computeCov(weights, weightsSum);
  }

  void
  computeMin(const Ref<const Matrix<T> >& data)
  {
    min_ = data.rowwise().minCoeff();
  }
  void
  computeMax(const Ref<const Matrix<T> >& data)
  {
    max_ = data.rowwise().maxCoeff();
  }
  void
  computeLimits(const Ref<const Matrix<T> >& data)
  {
    computeMin(data);
    computeMax(data);
  }
};
}  // namespace EventEMin

#endif  // EVENT_EMIN_DATA_STATS_H
