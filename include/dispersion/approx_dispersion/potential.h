#ifndef APPROX_POTENTIAL_H
#define APPROX_POTENTIAL_H

#include "dispersion/approx_dispersion/dispersion.h"

namespace event_model
{
template <typename M>
class ApproxPotential;

template <typename M>
struct DispersionTraits<ApproxDispersion<ApproxPotential<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template <typename M>
class ApproxPotential : public ApproxDispersion<ApproxPotential<M> >
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
  ApproxPotential(const T& dimScaleMax, const T& offset = T(5.0),
                  const T& lambda = T(0.0), const int ksize = 3)
      : ApproxDispersion<ApproxPotential<Model> >(dimScaleMax, offset, lambda,
                                                  ksize)
  {
  }

  template <typename U>
  U
  score(const Matrix<U>& c, const Convolution<U, NDims, T>& conv) const
  {
    return -this->interpolate(c, conv) / this->nPoints();
  }
};
}  // namespace event_model

#endif  // APPROX_POTENTIAL_H
