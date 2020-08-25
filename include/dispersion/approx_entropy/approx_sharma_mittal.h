#ifndef APPROX_SHARMA_MITTAL_H
#define APPROX_SHARMA_MITTAL_H

#include "dispersion/approx_entropy/approx_entropy.h"

namespace event_model
{

namespace dispersion
{

template<typename M>
class ApproxSharmaMittal;

template<typename M>
struct DispersionTraits<ApproxEntropy<ApproxSharmaMittal<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template<typename M>
class ApproxSharmaMittal : public ApproxEntropy<ApproxSharmaMittal<M> >
{

public:

  typedef M Model;
  typedef typename Model::T T;

private:

  const T gamma_;

protected:

  const T alpha_, beta_;

public:

  ApproxSharmaMittal(
    const array<Eigen::Index, Model::ndims()>& dim,
    const T& offset=T(5.0),
    const T& lambda=T(0.0),
    const int ksize=3,
    const T& alpha=T(2.0),
    const T& beta=T(0.5)):
    ApproxEntropy<ApproxSharmaMittal<Model> >(dim, offset, lambda, ksize),
    gamma_((1.0-beta)/(1.0-alpha)),
    alpha_(alpha), beta_(beta)
  {
    assert(0.0<alpha_ && alpha_!=1.0 && beta_!=1.0);
    this->kernel_=this->kernel_.pow(alpha_);
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

  template<typename U>
  U
  fFinal(
    const mtx<U>& c,
    const convolution::Convolution<U, Model::ndims(), T>& conv) const
  {
    return (pow((this->interpolate(c, conv)+T(1.0e-32))/this->npoints(), gamma())-T(1.0))/(T(1.0)-beta());
  }

protected:

private:

};

}

}

#endif // APPROX_SHARMA_MITTAL_H
