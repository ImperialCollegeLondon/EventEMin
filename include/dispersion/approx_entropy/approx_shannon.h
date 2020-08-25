#ifndef APPROX_SHANNON_H
#define APPROX_SHANNON_H

#include "dispersion/approx_entropy/approx_entropy.h"

namespace event_model
{

namespace dispersion
{

template<typename M>
class ApproxShannon;

template<typename M>
struct DispersionTraits<ApproxEntropy<ApproxShannon<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template<typename M>
class ApproxShannon : public ApproxEntropy<ApproxShannon<M> >
{

public:

  typedef M Model;
  typedef typename Model::T T;

private:

protected:

public:

  ApproxShannon(
    const array<Eigen::Index, Model::ndims()>& dim,
    const T& offset=T(5.0),
    const T& lambda=T(0.0),
    const int ksize=3):
    ApproxEntropy<ApproxShannon<Model> >(dim, offset, lambda, ksize)
  {
    this->kernel_=this->kernel_*this->kernel_.log();
  }

  template<typename U>
  U
  fFinal(
    const mtx<U>& c,
    const convolution::Convolution<U, Model::ndims(), T>& conv) const
  {
    return this->interpolate(c, conv)/this->npoints();
  }

protected:

private:

};

}

}

#endif // APPROX_SHANNON_H
