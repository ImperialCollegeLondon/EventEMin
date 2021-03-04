#ifndef TRANSLATION_NORMAL_H
#define TRANSLATION_NORMAL_H

#include "model/model.h"

namespace event_model
{

namespace model
{

namespace translation_normal
{

constexpr int nv=3, nvars=nv;
constexpr int ndims=2;
constexpr int nmtx=3;

template<typename Scalar>
struct TranslationNormal
{

  enum
  {
    InputsAtCompileTime=translation_normal::nvars,
    ValuesAtCompileTime=translation_normal::ndims
  };

  typedef Scalar T;
  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;

  TranslationNormal(void)
  {}

  static constexpr int
  nvars(void)
  {
    return translation_normal::nvars;
  }
  static constexpr int
  ndims(void)
  {
    return translation_normal::ndims;
  }
  static constexpr int
  nmtx(void)
  {
    return translation_normal::nmtx;
  }

  template<typename U>
  void
  operator()(
    const U* const vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    const Eigen::Map<const vec<U, translation_normal::nv> > v(vars);
    const vec<U, translation_normal::nv> vt(t*v);

    m.setIdentity();
    m.template rightCols<1>()-=vt;
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

    const Eigen::Map<const vec<V, translation_normal::ndims> > cMap(c);
    Eigen::Map<vec<U, translation_normal::ndims> > cmMap(cm);
    vec<V, translation_normal::ndims+1> ch;
    vec<U, translation_normal::ndims+1> chm;
    ch<<cMap, V(1.0);
    chm.noalias()=tMtx.inverse()*ch;
    cmMap=chm.template head<translation_normal::ndims>()/chm(translation_normal::ndims);
  }

};

}

template<typename T>
using TranslationNormal=Model<translation_normal::TranslationNormal<T> >;

}

}

#endif // TRANSLATION_NORMAL_H
