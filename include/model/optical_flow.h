#ifndef OPTICAL_FLOW_H
#define OPTICAL_FLOW_H

#include "model/model.h"

namespace event_model
{

namespace model
{

namespace optical_flow
{

constexpr int nv=2, nvars=nv;
constexpr int ndims=2;
constexpr int nmtx=3;

template<typename Scalar>
struct OpticalFlow
{

  enum
  {
    InputsAtCompileTime=optical_flow::nvars,
    ValuesAtCompileTime=optical_flow::ndims
  };

  typedef Scalar T;
  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;

  OpticalFlow(void)
  {}

  static constexpr int
  nvars(void)
  {
    return optical_flow::nvars;
  }
  static constexpr int
  ndims(void)
  {
    return optical_flow::ndims;
  }
  static constexpr int
  nmtx(void)
  {
    return optical_flow::nmtx;
  }

  template<typename U>
  void
  operator()(
    const U* const vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    const Eigen::Map<const vec<U, optical_flow::nv> > v(vars);

    m.setIdentity();
    m.template topRightCorner<optical_flow::nv, 1>()=t*v;
  }

  template<typename U, typename V>
  void
  operator()(
    const U* const vars,
    const V* const c,
    const T& t,
    U* cm) const
  {
    mtx<U, nmtx(), nmtx()> tMtx;
    (*this)(vars, t, tMtx);

    const Eigen::Map<const vec<V, optical_flow::ndims> > cMap(c);
    Eigen::Map<vec<U, optical_flow::ndims> > cmMap(cm);
    vec<V, optical_flow::ndims+1> ch;
    vec<U, optical_flow::ndims+1> chm;
    ch<<cMap, V(1.0);
    chm.noalias()=tMtx.inverse()*ch;
    cmMap=chm.template head<optical_flow::ndims>()/chm(optical_flow::ndims);
  }

};

}

template<typename T>
using OpticalFlow=Model<optical_flow::OpticalFlow<T> >;

}

}

#endif // OPTICAL_FLOW_H
