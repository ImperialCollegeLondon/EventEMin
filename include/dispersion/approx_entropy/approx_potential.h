#ifndef APPROX_POTENTIAL_H
#define APPROX_POTENTIAL_H

#include "dispersion/approx_entropy/approx_entropy.h"

namespace event_model
{

namespace dispersion
{

template<typename M>
class ApproxPotential;

template<typename M>
struct DispersionTraits<ApproxEntropy<ApproxPotential<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template<typename M>
class ApproxPotential : public ApproxEntropy<ApproxPotential<M> >
{

public:

  typedef M Model;
  typedef typename Model::T T;

private:

protected:

public:

  ApproxPotential(
    const array<Eigen::Index, Model::ndims()>& dim,
    const T& offset=T(5.0),
    const T& lambda=T(0.0),
    const int ksize=3):
    ApproxEntropy<ApproxPotential<Model> >(dim, offset, lambda, ksize)
  {}

  template<typename U>
  U
  fFinal(
    const mtx<U>& c,
    const convolution::Convolution<U, Model::ndims(), T>& conv) const
  {
    return -this->interpolate(c, conv)/this->npoints();
  }

protected:

private:

};

}

}

#endif // APPROX_POTENTIAL_H
