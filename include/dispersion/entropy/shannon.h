#ifndef SHANNON_H
#define SHANNON_H

#include "dispersion/entropy/entropy.h"

namespace event_model
{

namespace dispersion
{

template<typename M>
class Shannon;

template<typename M>
struct DispersionTraits<Entropy<Shannon<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template<typename M>
class Shannon : public Entropy<Shannon<M> >
{

public:

  typedef M Model;
  typedef typename Model::T T;

private:

protected:

  const T tau_, logTau_, scale_;

public:

  Shannon(
    const T& cScale=T(100.0)):
    Entropy<Shannon<Model> >(cScale),
    tau_(std::pow(2.0*M_PI, 0.5*this->ndims())),
    logTau_(std::log(tau_)),
    scale_(2.0/tau_)
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
    vec<U> cDiffPow, cDiffPowExp;
    pointsPowExp(cDiff, cDiffPow, cDiffPowExp);
    return (cDiffPowExp.array()*(cDiffPow.array()-logTau_)).sum();
  }

  template<typename U>
  U
  fFinal(
    const U& f) const
  {
    return scale()*(f/this->npoints());
  }

protected:

  template<typename U>
  void
  pointsPowExp(
    const mtx<U>& p,
    vec<U>& pPow,
    vec<U>& pPowExp) const
  {
    pPow=p.array().square().colwise().sum().transpose();
    pPowExp=(T(-0.5)*pPow.array()).exp();
  }

private:

};

}

}

#endif // SHANNON_H
