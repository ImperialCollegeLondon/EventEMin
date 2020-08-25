#ifndef TRANSLATION_3D_H
#define TRANSLATION_3D_H

#include "model/model.h"

namespace event_model
{

namespace model
{

namespace translation3d
{

constexpr int nv=3, nvars=nv;
constexpr int ndims=3;
constexpr int nmtx=4;

template<typename Scalar>
struct Translation3D
{

  enum
  {
    InputsAtCompileTime=translation3d::nvars,
    ValuesAtCompileTime=translation3d::ndims
  };

  typedef Scalar T;
  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;

  Translation3D(void)
  {}

  static constexpr int
  nvars(void)
  {
    return translation3d::nvars;
  }
  static constexpr int
  ndims(void)
  {
    return translation3d::ndims;
  }
  static constexpr int
  nmtx(void)
  {
    return translation3d::nmtx;
  }

  template<typename U>
  void
  operator()(
    const U* const vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    const Eigen::Map<const vec<U, translation3d::nvars> > varsMap(vars);
    const Eigen::Map<const vec<U, translation3d::nv> > v(varsMap.template segment<translation3d::nv>(0).data());

    mtx<U, nmtx(), nmtx()> skewMtx(mtx<U, nmtx(), nmtx()>::Zero());
    skewMtx.template rightCols<1>().template head<nv>()=v;
    m=mtx<U, nmtx(), nmtx()>::Identity()+t*skewMtx;
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

    const Eigen::Map<const vec<V, translation3d::ndims> > cMap(c);
    Eigen::Map<vec<U, translation3d::ndims> > cmMap(cm);
    vec<V, translation3d::ndims+1> ch;
    vec<U, translation3d::ndims+1> chm;
    ch<<cMap, V(1.0);
    chm.noalias()=tMtx.inverse()*ch;
    cmMap=chm.template head<translation3d::ndims>()/chm(translation3d::ndims);
  }

};

}

template<typename T>
using Translation3D=Model<translation3d::Translation3D<T> >;

}

}

#endif // TRANSLATION_3D_H
