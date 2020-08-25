#ifndef TSALLIS_H
#define TSALLIS_H

#include "dispersion/entropy/entropy.h"

namespace event_model
{

namespace dispersion
{

template<typename M>
class Tsallis;

template<typename M>
struct DispersionTraits<Entropy<Tsallis<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template<typename M>
class Tsallis : public Entropy<Tsallis<M> >
{

public:

  typedef M Model;
  typedef typename Model::T T;

private:

protected:

  const T alpha_, scale_;

public:

  Tsallis(
    const T& alpha=T(2.0),
    const T& cScale=T(100.0)):
    Entropy<Tsallis<Model> >(cScale),
    alpha_(alpha),
    scale_(2.0/std::pow(std::pow(2.0*M_PI, 0.5*this->ndims()), alpha_))
  {
    assert(0.0<alpha_ && alpha_!=1.0);
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
    return scale()*(f/this->npoints())/(T(1.0)-alpha());
  }

protected:

  template<typename U>
  void
  pointsPowExp(
    const mtx<U>& p,
    vec<U>& pPowExp) const
  {
    pPowExp=(T(-0.5*alpha())*p.array().square().colwise().sum()).exp();
  }

private:

};

}

}

#endif // TSALLIS_H
