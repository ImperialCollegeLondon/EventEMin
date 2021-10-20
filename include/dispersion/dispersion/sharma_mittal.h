#ifndef SHARMA_MITTAL_H
#define SHARMA_MITTAL_H

#include "dispersion/dispersion/dispersion.h"

namespace event_model
{
template <typename M>
class SharmaMittal;

template <typename M>
struct DispersionTraits<Dispersion<SharmaMittal<M> > >
{
  typedef M Model;
  typedef typename Model::T T;
};

template <typename M>
class SharmaMittal : public Dispersion<SharmaMittal<M> >
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
  const T alpha_, beta_, gamma_;
  const T scale_;

 public:
  SharmaMittal(const T& dimScale, const T& alpha = T(2.0),
               const T& beta = T(0.5))
      : Dispersion<SharmaMittal<Model> >(dimScale),
        alpha_(alpha),
        beta_(beta),
        gamma_((T(1.0) - beta_) / (T(1.0) - alpha_)),
        scale_(T(1.0) / (T(1.0) - beta_))
  {
    assert(T(0.0) < alpha_ && alpha_ != T(1.0));
    assert(beta_ != T(1.0));
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

  template <typename U>
  U
  partialScore(const Vector<U>& cDiffPow) const
  {
    return (alpha() * cDiffPow).array().exp().sum();
  }

  template <typename U>
  U
  score(const U& s) const
  {
    return scale() * pow(s / this->nPoints(), gamma());
  }
};
}  // namespace event_model

#endif  // SHARMA_MITTAL_H
