#ifndef DISPERSION_H
#define DISPERSION_H

#include <cassert>
#include <cmath>

#include "data_stats.h"
#include "gauss_kernel.h"
#include "types_def.h"

namespace event_model
{
template <typename Derived>
struct DispersionTraits;

template <typename Derived>
class DispersionBase
{
 public:
  typedef typename DispersionTraits<Derived>::Model Model;
  typedef typename DispersionTraits<Derived>::T T;

  enum
  {
    NVars = Model::NVars,
    NDims = Model::NDims,
    InputsAtCompileTime = NVars,
    ValuesAtCompileTime = 1
  };

  typedef Vector<T, InputsAtCompileTime> InputType;
  typedef Vector<T, ValuesAtCompileTime> ValueType;

 private:
  int nPoints_;

  DataStats<T> cStats_;
  Vector<T, NDims> cLimDiff_;

 protected:
  const Model model_;
  Map<const Matrix<T> > c_;
  Map<const Vector<T> > ts_;
  Map<const Vector<int> > polarity_;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  DispersionBase(void)
      : nPoints_(0), c_(nullptr, 0, 0), ts_(nullptr, 0), polarity_(nullptr, 0)
  {
  }

  int
  nPoints(void) const
  {
    return nPoints_;
  }
  const Ref<const Matrix<T> >
  c(void) const
  {
    return c_;
  }
  const Ref<const Vector<T> >
  c(const int k) const
  {
    assert(0 <= k && k < nPoints());
    return c_.col(k);
  }
  T
  cLimDiff(const int d) const
  {
    assert(0 <= d && d < NDims);
    return cLimDiff_(d);
  }
  const Vector<T, NDims>&
  cLimDiff(void) const
  {
    return cLimDiff_;
  }
  const Ref<const Vector<T> >
  ts(void) const
  {
    return ts_;
  }
  T
  ts(const int k) const
  {
    assert(0 <= k && k < nPoints());
    return ts_(k);
  }
  T
  tsRef(void) const
  {
    return ts_(0);
  }
  T
  tsMiddle(void) const
  {
    return ts_(nPoints() >> 1);
  }
  T
  tsEnd(void) const
  {
    return ts_(nPoints() - 1);
  }
  const Ref<const Vector<int> >
  polarity(void) const
  {
    return polarity_;
  }
  int
  polarity(const int k) const
  {
    assert(0 <= k && k < nPoints());
    return polarity_(k);
  }

  void
  assignPoints(const Ref<const Matrix<T> >& c, const Ref<const Vector<T> >& ts,
               const Ref<const Vector<int> >& polarity)
  {
    assert(c.rows() == NDims);
    nPoints_ = c.cols();
    assert(ts.size() == nPoints());
    assert(polarity.size() == nPoints());

    new (&c_)
        Map<const Matrix<T> >(c.data(), static_cast<int>(NDims), nPoints());
    new (&ts_) Map<const Vector<T> >(ts.data(), nPoints());
    new (&polarity_) Map<const Vector<int> >(polarity.data(), nPoints());

    cStats_.computeLimits(c_);
    cLimDiff_ = (cStats_.max() - cStats_.min()).array() + T(1.0e-8);

    this->underlying().computeDimScale();
  }

  template <typename U>
  void
  operator()(const Vector<U, NVars>& vars, Vector<U, 1>* f) const
  {
    Matrix<U> cm(static_cast<int>(NDims), nPoints()),
        cmScaled(static_cast<int>(NDims), nPoints());

    modelPoints(vars, cm);
    scalePoints<U>(cm, cmScaled);

    (*f)(0) = this->underlying().compute(cmScaled);
  }

 protected:
  template <typename U>
  void
  modelPoints(const Vector<U, NVars>& vars, Matrix<U>& cm) const
  {
    model_(vars, c(), ts(), cm);
  }

  template <typename U>
  void
  scalePoints(const Ref<const Matrix<U> >& c, Ref<Matrix<U> > cScaled) const
  {
    for (int d = 0; d < NDims; ++d)
    {
      cScaled.row(d) =
          ((this->underlying().dimScale(d) - this->underlying().offset()) /
           cLimDiff_(d)) *
              (c.row(d).array() - cStats_.min()(d)) +
          this->underlying().halfOffset();
    }
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
}  // namespace event_model

#endif  // DISPERSION_H
