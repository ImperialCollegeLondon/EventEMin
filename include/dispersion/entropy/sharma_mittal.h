#ifndef SHARMA_MITTAL_H
#define SHARMA_MITTAL_H

#include "dispersion/entropy/entropy.h"

namespace event_model
{

namespace dispersion
{

template<typename M>
class SharmaMittal;

template<typename M>
struct DispersionTraits<Entropy<SharmaMittal<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template<typename M>
class SharmaMittal : public Entropy<SharmaMittal<M> >
{

public:

  typedef M Model;
  typedef typename Model::T T;

private:

protected:

  const T alpha_, beta_, gamma_, scale_;

public:

  SharmaMittal(
    const T& alpha=T(2.0), const T& beta=T(0.5),
    const T& cScale=T(100.0)):
    Entropy<SharmaMittal<Model> >(cScale),
    alpha_(alpha), beta_(beta), gamma_((1.0-beta_)/(1.0-alpha_)),
    scale_(2.0/std::pow(std::pow(2.0*M_PI, 0.5*this->ndims()), alpha_))
  {
    assert(0.0<alpha_ && alpha_!=1.0 && beta_!=1.0);
  }

  T
  alpha(void) const
  {
    return alpha_;
  }
  T
  beta(void) const
  {
    return beta_;
  }
  T
  gamma(void) const
  {
    return gamma_;
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
    return (pow(scale()*((f+T(1.0e-32))/this->npoints()), gamma())-T(1.0))/(T(1.0)-beta());
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

#endif // SHARMA_MITTAL_H
