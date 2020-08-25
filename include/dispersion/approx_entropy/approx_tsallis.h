#ifndef APPROX_TSALLIS_H
#define APPROX_TSALLIS_H

#include "dispersion/approx_entropy/approx_entropy.h"

namespace event_model
{

namespace dispersion
{

template<typename M>
class ApproxTsallis;

template<typename M>
struct DispersionTraits<ApproxEntropy<ApproxTsallis<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template<typename M>
class ApproxTsallis : public ApproxEntropy<ApproxTsallis<M> >
{

public:

  typedef M Model;
  typedef typename Model::T T;

private:

protected:

  const T alpha_;

public:

  ApproxTsallis(
    const array<Eigen::Index, Model::ndims()>& dim,
    const T& offset=T(5.0),
    const T& lambda=T(0.0),
    const int ksize=3,
    const T& alpha=T(2.0)):
    ApproxEntropy<ApproxTsallis<Model> >(dim, offset, lambda, ksize),
    alpha_(alpha)
  {
    assert(0.0<alpha_ && alpha_!=1.0);
    this->kernel_=this->kernel_.pow(alpha_);
  }

  T
  alpha(void) const
  {
    return alpha_;
  }

  template<typename U>
  U
  fFinal(
    const mtx<U>& c,
    const convolution::Convolution<U, Model::ndims(), T>& conv) const
  {
    return this->interpolate(c, conv)/(this->npoints()*(T(1.0)-alpha()));
  }

protected:

private:

};

}

}

#endif // APPROX_TSALLIS_H
