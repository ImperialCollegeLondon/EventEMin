#ifndef EVENT_EMIN_SHANNON_H
#define EVENT_EMIN_SHANNON_H

#include "EventEMin/dispersion/dispersion/dispersion.h"

namespace EventEMin
{
namespace batch
{
namespace exact
{
template <typename M>
class Shannon;
}

template <typename M>
struct DispersionTraits<exact::Dispersion<exact::Shannon<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

namespace exact
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

 private:
  const T logDen_;

 public:
  Shannon(const T& dimScale)
      : Dispersion<Shannon<Model> >(dimScale),
        logDen_(std::log(std::pow(T(2.0 * M_PI), T(0.5) * NDims)))
  {
  }
  template <typename U>
  U
  partialScore(const Vector<U>& cDiffPow) const
  {
    return (cDiffPow.array().exp() * (cDiffPow.array() - logDen_)).sum();
  }

  template <typename U>
  U
  score(const U& s) const
  {
    return s / this->nPoints();
  }
};
}  // namespace exact
}  // namespace batch

template <typename M>
using Shannon = batch::exact::Shannon<M>;
}  // namespace EventEMin

#endif  // EVENT_EMIN_SHANNON_H
