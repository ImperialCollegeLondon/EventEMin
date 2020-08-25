#ifndef SIMILARITY_H
#define SIMILARITY_H

#include "model/model.h"

namespace event_model
{

namespace model
{

namespace similarity
{

constexpr int nw=1, ns=1, nv=2, nvars=nw+ns+nv;
constexpr int ndims=2;
constexpr int nmtx=3;

template<typename Scalar>
struct Similarity
{

  enum
  {
    InputsAtCompileTime=similarity::nvars,
    ValuesAtCompileTime=similarity::ndims
  };

  typedef Scalar T;
  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;

  Similarity(void)
  {}

  static constexpr int
  nvars(void)
  {
    return similarity::nvars;
  }
  static constexpr int
  ndims(void)
  {
    return similarity::ndims;
  }
  static constexpr int
  nmtx(void)
  {
    return similarity::nmtx;
  }

  template<typename U>
  void
  operator()(
    const U* const vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    const Eigen::Map<const vec<U, similarity::nvars> > varsMap(vars);
    const vec<U, similarity::nvars> varst(t*varsMap);
    const Eigen::Map<const vec<U, similarity::nw> > wt(varst.template segment<similarity::nw>(0).data());
    const Eigen::Map<const vec<U, similarity::ns> > st(varst.template segment<similarity::ns>(similarity::nw).data());
    const Eigen::Map<const vec<U, similarity::nv> > vt(varst.template segment<similarity::nv>(similarity::nw+similarity::ns).data());

    const U costheta=cos(wt(0)), sintheta=sin(wt(0));

    mtx<U, 2, 2> rthetaMtx;
    rthetaMtx<<costheta, -sintheta,
      sintheta, costheta;

    const mtx<U, 2, 2> sMtx((st(0)+U(1.0))*rthetaMtx);
    m<<sMtx(0, 0), sMtx(0, 1), vt(0),
      sMtx(1, 0), sMtx(1, 1), vt(1),
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

    const Eigen::Map<const vec<V, similarity::ndims> > cMap(c);
    Eigen::Map<vec<U, similarity::ndims> > cmMap(cm);
    vec<V, similarity::ndims+1> ch;
    vec<U, similarity::ndims+1> chm;
    ch<<cMap, V(1.0);
    chm.noalias()=tMtx.inverse()*ch;
    cmMap=chm.template head<similarity::ndims>()/chm(similarity::ndims);
  }

};

}

template<typename T>
using Similarity=Model<similarity::Similarity<T> >;

}

}

#endif // SIMILARITY_H
