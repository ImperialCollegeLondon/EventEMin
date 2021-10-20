#ifndef APPROX_SHANNON_H
#define APPROX_SHANNON_H

#include "dispersion/approx_dispersion/dispersion.h"

namespace event_model
{
template <typename M>
class ApproxShannon;

template <typename M>
struct DispersionTraits<ApproxDispersion<ApproxShannon<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template <typename M>
class ApproxShannon : public ApproxDispersion<ApproxShannon<M> >
{
 public:
  typedef M Model;
  typedef typename Model::T T;

  enum
  {
    NVars = Model::NVars,
    NDims = Model::NDims
  };

 public:
  ApproxShannon(const T& dimScaleMax, const T& offset = T(5.0),
                const T& lambda = T(0.0), const int ksize = 3)
      : ApproxDispersion<ApproxShannon<Model> >(dimScaleMax, offset, lambda,
                                                ksize)
  {
    this->kernel_ = this->kernel_ * this->kernel_.log();
  }

  template <typename U>
  U
  score(const Matrix<U>& c, const Convolution<U, NDims, T>& conv) const
  {
    return this->interpolate(c, conv) / this->nPoints();
  }
};
}  // namespace event_model

#endif  // APPROX_SHANNON_H
