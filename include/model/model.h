#ifndef MODEL_H
#define MODEL_H

#include <cassert>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>

#include "types_def.h"

namespace event_model
{

template<typename F>
class Model
{

public:

  typedef F Func;
  typedef typename Func::T T;

  enum
  {
    InputsAtCompileTime=Func::nvars(),
    ValuesAtCompileTime=Func::ndims()
  };

  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;
  
private:

protected:

  const Func func_;

public:

  static constexpr int
  nvars(void)
  {
    return Func::nvars();
  }
  static constexpr int
  ndims(void)
  {
    return Func::ndims();
  }
  static constexpr int
  nmtx(void)
  {
    return Func::nmtx();
  }

  template<typename U>
  void
  operator()(
    const U* const vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    func_(vars, t, m);
  }

  template<typename U>
  void
  operator()(
    const vec<U, nvars()>& vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    (*this)(vars.data(), t, m);
  }

  template<typename U, typename V>
  void
  operator()(
    const U* const vars,
    const V* const c,
    const T& t,
    U* cm) const
  {
    func_(vars, c, t, cm);
  }

  template<typename U>
  void
  operator()(
    const vec<U, nvars()>& vars,
    vec<U, ndims()>* cm,
    const Eigen::Ref<const vec<T> >& c,
    const T& t) const
  {
    assert(c.size()==ndims());
    
    (*this)(vars.data(), c.data(), t, cm->data());
  }

  template<typename U>
  void
  operator()(
    const vec<U, nvars()>& vars,
    const Eigen::Ref<const mtx<T> >& c,
    const Eigen::Ref<const vec<T> >& ts,
    const T& tsref,
    mtx<U>& cm) const
  {
    assert(c.rows()==ndims() && cm.rows()==ndims());
    assert(c.cols()==ts.size() && c.cols()==cm.cols());

    const vec<T> t(ts.array()-tsref);
    vec<U, ndims()> vcm;
    for(int k=0; k<c.cols(); ++k)
    {
      (*this)(vars, &vcm, c.col(k), t(k));
      cm.col(k)=vcm;
    }
  }

  template<typename U>
  void
  operator()(
    const vec<U, nvars()>& vars,
    const Eigen::Ref<const mtx<T> >& c,
    const Eigen::Ref<const vec<T> >& ts,
    mtx<U>& cm) const
  {
    (*this)(vars, c, ts, ts(0), cm);
  }

protected:

private:

};

}

#endif // MODEL_H
