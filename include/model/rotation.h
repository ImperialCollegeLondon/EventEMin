#ifndef ROTATION_H
#define ROTATION_H

#include "model/model.h"

namespace event_model
{

namespace model
{

namespace rotation
{

constexpr int nw=3, nvars=nw;
constexpr int ndims=2;
constexpr int nmtx=3;

template<typename Scalar>
struct Rotation
{

  enum
  {
    InputsAtCompileTime=rotation::nvars,
    ValuesAtCompileTime=rotation::ndims
  };

  typedef Scalar T;
  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;

  Rotation(void)
  {}

  static constexpr int
  nvars(void)
  {
    return rotation::nvars;
  }
  static constexpr int
  ndims(void)
  {
    return rotation::ndims;
  }
  static constexpr int
  nmtx(void)
  {
    return rotation::nmtx;
  }

  template<typename U>
  void
  operator()(
    const U* const vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    const Eigen::Map<const vec<U, rotation::nw> > w(vars);

    m.setIdentity();
    const U wNorm=w.norm();
    if(U(0.0)<wNorm)
    {
      const U theta=t*wNorm;
      mtx<U, nmtx(), nmtx()> skewMtx;
      skewMtx<<U(0.0), -w(2), w(1),
        w(2), U(0.0), -w(0),
        -w(1), w(0), U(0.0);
      skewMtx/=wNorm;
      const U ctheta=U(1.0)-cos(theta);
      const U stheta=sin(theta);
      m+=stheta*skewMtx+ctheta*skewMtx*skewMtx;
    }
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

    const Eigen::Map<const vec<V, rotation::ndims> > cMap(c);
    Eigen::Map<vec<U, rotation::ndims> > cmMap(cm);
    vec<V, rotation::ndims+1> ch;
    vec<U, rotation::ndims+1> chm;
    ch<<cMap, V(1.0);
    chm.noalias()=tMtx.transpose()*ch;
    cmMap=chm.template head<rotation::ndims>()/chm(rotation::ndims);
  }

};

}

template<typename T>
using Rotation=Model<rotation::Rotation<T> >;

}

}

#endif // ROTATION_H
