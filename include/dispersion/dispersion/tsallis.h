#ifndef TSALLIS_H
#define TSALLIS_H

#include "dispersion/dispersion/dispersion.h"

namespace event_model
{
template <typename M>
class Tsallis;

template <typename M>
struct DispersionTraits<Dispersion<Tsallis<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template <typename M>
class Tsallis : public Dispersion<Tsallis<M> >
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
  Tsallis(const T& dimScale, const T& alpha = T(2.0))
      : Dispersion<Tsallis<Model> >(dimScale),
        alpha_(alpha),
        scale_(T(1.0) / (T(1.0) - alpha_))
  {
    assert(T(0.0) < alpha_ && alpha_ != T(1.0));
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
  partialScore(const Vector<U>& cDiffPow) const
  {
    return (alpha() * cDiffPow).array().exp().sum();
  }

  template <typename U>
  U
  score(const U& s) const
  {
    return scale() * s / this->nPoints();
  }
};
}  // namespace event_model

#endif  // TSALLIS_H
