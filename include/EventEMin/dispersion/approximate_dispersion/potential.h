#ifndef EVENT_EMIN_APPROXIMATE_POTENTIAL_H
#define EVENT_EMIN_APPROXIMATE_POTENTIAL_H

#include "EventEMin/dispersion/approximate_dispersion/dispersion.h"

namespace EventEMin
{
namespace batch
{
namespace approximate
{
template <typename M>
class Potential;
}

template <typename M>
struct DispersionTraits<approximate::Dispersion<approximate::Potential<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

namespace approximate
{
template <typename M>
class Potential : public Dispersion<Potential<M> >
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
  Potential(const T& dimScaleMax, const T& offset = T(5.0),
            const T& lambda = T(0.0), const int ksize = 3)
      : Dispersion<Potential<Model> >(dimScaleMax, offset, lambda, ksize)
  {
  }

  template <typename U>
  U
  score(const Matrix<U>& c, const Convolution<U, NDims, T>& conv) const
  {
    return -this->interpolate(c, conv) / this->nPoints();
  }
};
}  // namespace approximate
}  // namespace batch

template <typename M>
using ApproximatePotential = batch::approximate::Potential<M>;
}  // namespace EventEMin

#endif  // EVENT_EMIN_APPROXIMATE_POTENTIAL_H
