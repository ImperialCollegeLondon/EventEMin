#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "dispersion/entropy/entropy.h"

namespace event_model
{

namespace dispersion
{

template<typename M>
class Potential;

template<typename M>
struct DispersionTraits<Entropy<Potential<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template<typename M>
class Potential : public Entropy<Potential<M> >
{

public:

  typedef M Model;
  typedef typename Model::T T;

private:

protected:

  const T scale_;

public:

  Potential(
    const T& cScale=T(100.0)):
    Entropy<Potential<Model> >(cScale),
    scale_(2.0/(std::pow(2.0*M_PI, 0.5*this->ndims())))
  {}

  T
  scale(void) const
  {
    return scale_;
  }

  template<typename U>
  U
  fCumulative(
    const mtx<U>& cDiff) const
  {
    vec<U> cDiffPowExp;
    pointsPowExp(cDiff, cDiffPowExp);
    return cDiffPowExp.sum();
  }

  template<typename U>
  U
  fFinal(
    const U& f) const
  {
    return -scale()*(f/this->npoints());
  }

protected:

  template<typename U>
  void
  pointsPowExp(
    const mtx<U>& p,
    vec<U>& pPowExp) const
  {
    pPowExp=(T(-0.5)*p.array().square().colwise().sum()).exp();
  }

private:

};

}

}

#endif // POTENTIAL_H
