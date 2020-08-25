#ifndef TRANSLATION_H
#define TRANSLATION_H

#include "model/model.h"

namespace event_model
{

namespace model
{

namespace translation
{

constexpr int nv=3, nnsph=2, nn=3, nvars=nv+nnsph;
constexpr int ndims=2;
constexpr int nmtx=3;

template<typename Scalar>
struct Translation
{

  enum
  {
    InputsAtCompileTime=translation::nvars,
    ValuesAtCompileTime=translation::ndims
  };

  typedef Scalar T;
  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;

  Translation(void)
  {}

  static constexpr int
  nvars(void)
  {
    return translation::nvars;
  }
  static constexpr int
  ndims(void)
  {
    return translation::ndims;
  }
  static constexpr int
  nmtx(void)
  {
    return translation::nmtx;
  }

  template<typename U>
  void
  operator()(
    const U* const vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    const Eigen::Map<const vec<U, translation::nvars> > varsMap(vars);
    const vec<U, translation::nv> varst(t*varsMap.template head<translation::nv>());
    const Eigen::Map<const vec<U, translation::nv> > vt(varst.template segment<translation::nv>(0).data());
    const Eigen::Map<const vec<U, translation::nnsph> > nsph(varsMap.template segment<translation::nnsph>(translation::nv).data());

    const U costheta=cos(nsph(0)), sintheta=sin(nsph(0));
    const U cosphi=cos(nsph(1)), sinphi=sin(nsph(1));
    vec<U, translation::nn> n;
    n<<costheta*sinphi, sintheta*sinphi, cosphi;

    m=mtx<U, nmtx(), nmtx()>::Identity()-vt*n.transpose();
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

    const Eigen::Map<const vec<V, translation::ndims> > cMap(c);
    Eigen::Map<vec<U, translation::ndims> > cmMap(cm);
    vec<V, translation::ndims+1> ch;
    vec<U, translation::ndims+1> chm;
    ch<<cMap, V(1.0);
    chm.noalias()=tMtx.inverse()*ch;
    cmMap=chm.template head<translation::ndims>()/chm(translation::ndims);
  }

};

}

template<typename T>
using Translation=Model<translation::Translation<T> >;

}

}

#endif // TRANSLATION_H
