#ifndef EVENT_EMIN_APPROXIMATE_RENYI_H
#define EVENT_EMIN_APPROXIMATE_RENYI_H

#include "EventEMin/dispersion/approximate_dispersion/dispersion.h"

namespace EventEMin
{
namespace batch
{
namespace approximate
{
template <typename M>
class Renyi;
}

template <typename M>
struct DispersionTraits<approximate::Dispersion<approximate::Renyi<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

namespace approximate
{
template <typename M>
class Renyi : public Dispersion<Renyi<M> >
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
  Renyi(const T& dimScaleMax, const T& offset = T(5.0),
        const T& lambda = T(0.0), const int ksize = 3, const T& alpha = T(2.0))
      : Dispersion<Renyi<Model> >(dimScaleMax, offset, lambda, ksize),
        alpha_(alpha),
        scale_(T(1.0) / (T(1.0) - alpha_))
  {
    assert(T(0.0) < alpha_ && alpha_ != T(1.0));
    if constexpr (NDims == 2)
    {
      this->kernel_ = this->kernel_.array().pow(alpha_);
    }
    else
    {
      this->kernel_ = this->kernel_.pow(alpha_);
    }
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
    return scale() * log(this->interpolate(c, conv) / this->nPoints());
  }
};
}  // namespace approximate
}  // namespace batch

template <typename M>
using ApproximateRenyi = batch::approximate::Renyi<M>;
}  // namespace EventEMin

#endif  // EVENT_EMIN_APPROXIMATE_RENYI_H
