#ifndef EVENT_EMIN_APPROXIMATE_SHANNON_H
#define EVENT_EMIN_APPROXIMATE_SHANNON_H

#include "EventEMin/dispersion/approximate_dispersion/dispersion.h"

namespace EventEMin
{
namespace batch
{
namespace approximate
{
template <typename M>
class Shannon;
}

template <typename M>
struct DispersionTraits<approximate::Dispersion<approximate::Shannon<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

namespace approximate
{
template <typename M>
class Shannon : public Dispersion<Shannon<M> >
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
  Shannon(const T& dimScaleMax, const T& offset = T(5.0),
          const T& lambda = T(0.0), const int ksize = 3)
      : Dispersion<Shannon<Model> >(dimScaleMax, offset, lambda, ksize)
  {
    if constexpr (NDims == 2)
    {
      this->kernel_ = this->kernel_.array() * this->kernel_.array().log();
    }
    else
    {
      this->kernel_ = this->kernel_ * this->kernel_.log();
    }
  }

  template <typename U>
  U
  score(const Matrix<U>& c, const Convolution<U, NDims, T>& conv) const
  {
    return this->interpolate(c, conv) / this->nPoints();
  }
};
}  // namespace approximate
}  // namespace batch

template <typename M>
using ApproximateShannon = batch::approximate::Shannon<M>;
}  // namespace EventEMin

#endif  // EVENT_EMIN_APPROX_SHANNON_H
