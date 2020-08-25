#ifndef APPROX_RENYI_H
#define APPROX_RENYI_H

#include "dispersion/approx_entropy/approx_entropy.h"

namespace event_model
{

namespace dispersion
{

template<typename M>
class ApproxRenyi;

template<typename M>
struct DispersionTraits<ApproxEntropy<ApproxRenyi<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template<typename M>
class ApproxRenyi : public ApproxEntropy<ApproxRenyi<M> >
{

public:

  typedef M Model;
  typedef typename Model::T T;

private:

protected:

  const T alpha_;

public:

  ApproxRenyi(
    const array<Eigen::Index, Model::ndims()>& dim,
    const T& offset=T(5.0),
    const T& lambda=T(0.0),
    const int ksize=3,
    const T& alpha=T(2.0)):
    ApproxEntropy<ApproxRenyi<Model> >(dim, offset, lambda, ksize),
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
    return log((this->interpolate(c, conv)+T(1.0e-32))/this->npoints())/(T(1.0)-alpha());
  }

protected:

private:

};

}

}

#endif // APPROX_RENYI_H
