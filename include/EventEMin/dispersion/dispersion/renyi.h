#ifndef EVENT_EMIN_RENYI_H
#define EVENT_EMIN_RENYI_H

#include "EventEMin/dispersion/dispersion/dispersion.h"

namespace EventEMin
{
namespace batch
{
namespace exact
{
template <typename M>
class Renyi;
}

template <typename M>
struct DispersionTraits<exact::Dispersion<exact::Renyi<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

namespace exact
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
  Renyi(const T& dimScale, const T& alpha = T(2.0))
      : Dispersion<Renyi<Model> >(dimScale),
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
    return scale() * log(s / this->nPoints());
  }
};
}  // namespace exact
}  // namespace batch

template <typename M>
using Renyi = batch::exact::Renyi<M>;
}  // namespace EventEMin

#endif  // EVENT_EMIN_RENYI_H
