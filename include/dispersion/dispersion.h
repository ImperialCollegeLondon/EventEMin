#ifndef DISPERSION_H
#define DISPERSION_H

#include <cassert>
#include <cmath>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>

#include "data_stats.h"
#include "gauss_kernel.h"
#include "types_def.h"

namespace event_model
{

namespace dispersion
{

template<typename Derived>
struct DispersionTraits;

template<typename Derived>
class Dispersion
{

public:

  typedef typename DispersionTraits<Derived>::Model Model;
  typedef typename DispersionTraits<Derived>::T T;

  enum
  {
    InputsAtCompileTime=Model::nvars(),
    ValuesAtCompileTime=1
  };

  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;

private:

  int npoints_;

  T climDiffMax_;
  DataStats<T> cStats_;

protected:

  const Model model_;
  Eigen::Map<const mtx<T> > c_;
  Eigen::Map<const vec<T> > ts_;
  Eigen::Map<const vec<int> > polarity_;

public:

  Dispersion(void):
    npoints_(0),
    c_(nullptr, 0, 0), ts_(nullptr, 0), polarity_(nullptr, 0)
  {}

  static constexpr int
  nvars(void)
  {
    return Model::nvars();
  }
  static constexpr int
  ndims(void)
  {
    return Model::ndims();
  }
  int
  npoints(void) const
  {
    return npoints_;
  }
  const Eigen::Ref<const mtx<T> >
  c(void) const
  {
    return c_;
  }
  const Eigen::Ref<const vec<T> >
  climMin(void) const
  {
    return cStats_.minimum();
  }
  T
  climMin(
    const int d) const
  {
    return cStats_.minimum()(d);
  }
  const Eigen::Ref<const vec<T> >
  climMax(void) const
  {
    return cStats_.maximum();
  }
  T
  climMax(
    const int d) const
  {
    return cStats_.maximum()(d);
  }
  T
  climDiffMax(void) const
  {
    return climDiffMax_;
  }
  const Eigen::Ref<const vec<T> >
  ts(void) const
  {
    return ts_;
  }
  T
  ts(
    const int k) const
  {
    return ts_(k);
  }
  T
  tsref(void) const
  {
    return ts_(0);
  }
  T
  tsend(void) const
  {
    return ts_.template tail<1>()(0);
  }
  const Eigen::Ref<const vec<int> >
  polarity(void) const
  {
    return polarity_;
  }
  int
  polarity(
    const int k) const
  {
    return polarity_(k);
  }

  void
  assignPoints(
    const Eigen::Ref<const mtx<T> >& c,
    const Eigen::Ref<const vec<T> >& ts,
    const Eigen::Ref<const vec<int> >& polarity)
  {
    assert(ndims()==c.rows());
    npoints_=c.cols();
    assert(npoints()==ts.size() && npoints()==polarity.size());
    new (&c_) Eigen::Map<const mtx<T> >(c.data(), ndims(), npoints());
    new (&ts_) Eigen::Map<const vec<T> >(ts.data(), npoints());
    new (&polarity_) Eigen::Map<const vec<int> >(polarity.data(), npoints());

    cStats_.computeAll(c_);
    climDiffMax_=(climMax()-climMin()).maxCoeff()+T(1.0e-8);
  }

  template<typename U>
  void
  operator()(
    const vec<U, nvars()>& vars,
    vec<U, 1>* f) const
  {
    mtx<U> cm(ndims(), npoints()), cmScaled(ndims(), npoints());

    modelPoints(vars, cm);
    scalePoints(cm, cmScaled);

    (*f)(0)=this->underlying().compute(cmScaled);
  }

protected:

  template<typename U>
  void
  modelPoints(
    const vec<U, nvars()>& vars,
    mtx<U>& cm) const
  {
    model_(vars, c(), ts(), cm);
  }

  template<typename U>
  void
  scalePoints(
    const mtx<U>& c,
    mtx<U>& cScaled) const
  {
    for(int d=0; d<ndims(); ++d)
    {
      cScaled.row(d)=((this->underlying().dimScale(d)-this->underlying().offset())/climDiffMax())*(c.row(d).array()-climMin(d))+this->underlying().halfOffset();
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

}

}

#endif // DISPERSION_H
