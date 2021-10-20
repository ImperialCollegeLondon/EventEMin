#ifndef SHANNON_H
#define SHANNON_H

#include "dispersion/dispersion/dispersion.h"

namespace event_model
{
template <typename M>
class Shannon;

template <typename M>
struct DispersionTraits<Dispersion<Shannon<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

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
}  // namespace event_model

#endif  // SHANNON_H
