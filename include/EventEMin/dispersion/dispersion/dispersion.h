#ifndef EVENT_EMIN_DISPERSION_IMPL_H
#define EVENT_EMIN_DISPERSION_IMPL_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include "EventEMin/dispersion/dispersion.h"

namespace EventEMin
{
namespace batch
{
namespace exact
{
template <typename Derived>
class Dispersion : public DispersionBase<Dispersion<Derived> >
{
 public:
  typedef typename DispersionTraits<Dispersion<Derived> >::Model Model;
  typedef typename DispersionTraits<Dispersion<Derived> >::T T;

  enum
  {
    NVars = Model::NVars,
    NDims = Model::NDims
  };

 private:
  const T dimScaleMax_;
  Array<Index, NDims> dim_;

 public:
  Dispersion(const T& dimScaleMax)
      : DispersionBase<Dispersion<Derived> >(), dimScaleMax_(dimScaleMax)
  {
  }

  T
  dimScale(const int d) const
  {
    assert(0 <= d && d < NDims);
    return static_cast<T>(dim_[d]);
  }
  T
  halfOffset(void) const
  {
    return T(0.0);
  }
  T
  offset(void) const
  {
    return T(0.0);
  }

  template <typename U>
  U
  compute(const Matrix<U>& c) const
  {
    U s = U(0.0);
#ifdef _OPENMP
#pragma omp parallel shared(c, s)
#endif
    {
      Matrix<U> cDiff;
      Vector<U> cDiffPow;

#ifdef _OPENMP
#pragma omp declare reduction(sum:U : omp_out += omp_in) initializer(omp_priv = U(0.0))
#pragma omp for reduction(sum : s)
#endif
      for (int k = 0; k < this->nPoints(); ++k)
      {
        pointsDiff(k, c, cDiff);
        pointsPow(cDiff, cDiffPow);
        s += this->underlying().partialScore(cDiffPow);
      }
    }

    return this->underlying().score(s);
  }

  void
  computeDimScale(void)
  {
    const Vector<T, NDims> cLimDiff(this->c().rowwise().maxCoeff() -
                                    this->c().rowwise().minCoeff());
    const T limDiffMax = cLimDiff.maxCoeff();
    for (int d = 0; d < NDims; ++d)
    {
      dim_[d] = std::round(dimScaleMax_ * cLimDiff(d) / limDiffMax);
    }
  }

 protected:
  template <typename U>
  void
  pointsDiff(const int k, const Matrix<U>& p, Matrix<U>& pDiff) const
  {
    pDiff = p.colwise() - p.col(k);
  }

  template <typename U>
  void
  pointsPow(const Matrix<U>& p, Vector<U>& pPow) const
  {
    pPow = T(-0.5) * p.array().square().colwise().sum();
  }

 private:
  Derived&
  underlying(void)
  {
    return static_cast<Derived&>(*this);
  }
  const Derived&
  underlying(void) const
  {
    return static_cast<const Derived&>(*this);
  }
};
}  // namespace exact
}  // namespace batch
}  // namespace EventEMin

#endif  // EVENT_EMIN_DISPERSION_IMPL_H
