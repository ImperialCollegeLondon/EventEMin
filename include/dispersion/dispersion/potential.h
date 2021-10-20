#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "dispersion/dispersion/dispersion.h"

namespace event_model
{
template <typename M>
class Potential;

template <typename M>
struct DispersionTraits<Dispersion<Potential<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

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
}  // namespace event_model

#endif  // POTENTIAL_H
