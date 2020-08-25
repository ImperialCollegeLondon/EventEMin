#ifndef DATA_STATS_H
#define DATA_STATS_H

#include <cassert>
#include <Eigen/Core>

#include "types_def.h"

template<typename T>
class DataStats
{

private:

protected:

  vec<T> minimum_, maximum_;
  vec<T> mean_;
  mtx<T> centred_, covariance_;

public:

  const vec<T>&
  minimum(void) const
  {
    return minimum_;
  }
  const vec<T>&
  maximum(void) const
  {
    return maximum_;
  }
  const vec<T>&
  mean(void) const
  {
    return mean_;
  }
  const mtx<T>&
  centred(void) const
  {
    return centred_;
  }
  const mtx<T>&
  covariance(void) const
  {
    return covariance_;
  }

  void
  computeAll(
    const Eigen::Ref<const mtx<T> >& data)
  {
    computeLimits(data);
    computeMoments(data);
  }
  void
  computeAll(
    const Eigen::Ref<const mtx<T> >& data,
    const Eigen::Ref<const vec<T> >& weights)
  {
    assert(data.cols()==weights.size());
    computeLimits(data);
    computeMoments(data, weights, weights.sum()+T(1.0e-32));
  }

  void
  computeCentred(
    const Eigen::Ref<const mtx<T> >& data)
  {
    centred_.noalias()=data.colwise()-mean();
  }
  void
  computeCovariance(void)
  {
    assert(1<centred().cols());
    covariance_.noalias()=centred()*centred().transpose()/(centred().cols()-1);
  }
  void
  computeCovariance(
    const Eigen::Ref<const vec<T> >& weights,
    const T& weightsSum)
  {
    covariance_.noalias()=(centred().array().rowwise()*weights.transpose().array()).matrix()*centred().transpose()/weightsSum;
  }
  void
  computeMean(
    const Eigen::Ref<const mtx<T> >& data)
  {
    mean_.noalias()=data.rowwise().mean();
  }
  void
  computeMean(
    const Eigen::Ref<const mtx<T> >& data,
    const Eigen::Ref<const vec<T> >& weights,
    const T& weightsSum)
  {
    mean_.noalias()=(data.array().rowwise()*weights.transpose().array()).matrix().rowwise().sum()/weightsSum;
  }
  void
  computeMoments(
    const Eigen::Ref<const mtx<T> >& data)
  {
    computeMean(data);
    computeCentred(data);
    computeCovariance();
  }
  void
  computeMoments(
    const Eigen::Ref<const mtx<T> >& data,
    const Eigen::Ref<const vec<T> >& weights,
    const T& weightsSum)
  {
    computeMean(data, weights, weightsSum);
    computeCentred(data);
    computeCovariance(weights, weightsSum);
  }

  void
  computeMinimum(
    const Eigen::Ref<const mtx<T> >& data)
  {
    minimum_=data.rowwise().minCoeff();
  }
  void
  computeMaximum(
    const Eigen::Ref<const mtx<T> >& data)
  {
    maximum_=data.rowwise().maxCoeff();
  }
  void
  computeLimits(
    const Eigen::Ref<const mtx<T> >& data)
  {
    computeMinimum(data);
    computeMaximum(data);
  }

protected:

private:

};

#endif // DATA_STATS_H
