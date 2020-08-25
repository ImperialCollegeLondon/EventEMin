#ifndef ISOMETRY_H
#define ISOMETRY_H

#include "model/model.h"

namespace event_model
{

namespace model
{

namespace isometry
{

constexpr int nw=1, nv=2, nvars=nw+nv;
constexpr int ndims=2;
constexpr int nmtx=3;

template<typename Scalar>
struct Isometry
{

  enum
  {
    InputsAtCompileTime=isometry::nvars,
    ValuesAtCompileTime=isometry::ndims
  };

  typedef Scalar T;
  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;

  Isometry(void)
  {}

  static constexpr int
  nvars(void)
  {
    return isometry::nvars;
  }
  static constexpr int
  ndims(void)
  {
    return isometry::ndims;
  }
  static constexpr int
  nmtx(void)
  {
    return isometry::nmtx;
  }

  template<typename U>
  void
  operator()(
    const U* const vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    const Eigen::Map<const vec<U, isometry::nvars> > varsMap(vars);
    const vec<U, isometry::nvars> varst(t*varsMap);
    const Eigen::Map<const vec<U, isometry::nw> > wt(varst.template segment<isometry::nw>(0).data());
    const Eigen::Map<const vec<U, isometry::nv> > vt(varst.template segment<isometry::nv>(isometry::nw).data());

    const U costheta=cos(wt(0)), sintheta=sin(wt(0));

    mtx<U, 2, 2> rthetaMtx;
    rthetaMtx<<costheta, -sintheta,
      sintheta, costheta;

    m<<rthetaMtx(0, 0), rthetaMtx(0, 1), vt(0),
      rthetaMtx(1, 0), rthetaMtx(1, 1), vt(1),
      U(0.0), U(0.0), U(1.0);
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

    const Eigen::Map<const vec<V, isometry::ndims> > cMap(c);
    Eigen::Map<vec<U, isometry::ndims> > cmMap(cm);
    vec<V, isometry::ndims+1> ch;
    vec<U, isometry::ndims+1> chm;
    ch<<cMap, V(1.0);
    chm.noalias()=tMtx.inverse()*ch;
    cmMap=chm.template head<isometry::ndims>()/chm(isometry::ndims);
  }

};

}

template<typename T>
using Isometry=Model<isometry::Isometry<T> >;

}

}

#endif // ISOMETRY_H
