#ifndef APPROX_TSALLIS_H
#define APPROX_TSALLIS_H

#include "dispersion/approx_dispersion/dispersion.h"

namespace event_model
{
template <typename M>
class ApproxTsallis;

template <typename M>
struct DispersionTraits<ApproxDispersion<ApproxTsallis<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template <typename M>
class ApproxTsallis : public ApproxDispersion<ApproxTsallis<M> >
{
 public:
  typedef M Model;
  typedef typename Model::T T;

  enum
  {
    NVars = Model::NVars,
    NDims = Model::NDims
  };

 private:
  const T alpha_;
  const T scale_;

 public:
  ApproxTsallis(const T& dimScaleMax, const T& offset = T(5.0),
                const T& lambda = T(0.0), const int ksize = 3,
                const T& alpha = T(2.0))
      : ApproxDispersion<ApproxTsallis<Model> >(dimScaleMax, offset, lambda,
                                                ksize),
        alpha_(alpha),
        scale_(T(1.0) / (T(1.0) - alpha_))
  {
    assert(T(0.0) < alpha_ && alpha_ != T(1.0));
    this->kernel_ = this->kernel_.pow(alpha_);
  }

  T
  alpha(void) const
  {
    return alpha_;
  }
  T
  scale(void) const
  {
    return scale_;
  }

  template <typename U>
  U
  score(const Matrix<U>& c, const Convolution<U, NDims, T>& conv) const
  {
    return scale() * this->interpolate(c, conv) / this->nPoints();
  }
};
}  // namespace event_model

#endif  // APPROX_TSALLIS_H
