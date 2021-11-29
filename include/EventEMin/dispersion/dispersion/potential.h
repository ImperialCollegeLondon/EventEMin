#ifndef EVENT_EMIN_POTENTIAL_H
#define EVENT_EMIN_POTENTIAL_H

#include "EventEMin/dispersion/dispersion/dispersion.h"

namespace EventEMin
{
namespace batch
{
namespace exact
{
template <typename M>
class Potential;
}

template <typename M>
struct DispersionTraits<exact::Dispersion<exact::Potential<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

namespace exact
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
  Potential(const T& dimScale) : Dispersion<Potential<Model> >(dimScale) {}
  template <typename U>
  U
  partialScore(const Vector<U>& cDiffPow) const
  {
    return cDiffPow.array().exp().sum();
  }

  template <typename U>
  U
  score(const U& s) const
  {
    return -s / this->nPoints();
  }
};
}  // namespace exact
}  // namespace batch

template <typename M>
using Potential = batch::exact::Potential<M>;
}  // namespace EventEMin

#endif  // EVENT_EMIN_POTENTIAL_H
