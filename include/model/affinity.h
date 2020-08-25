#ifndef AFFINITY_H
#define AFFINITY_H

#include "model/model.h"

namespace event_model
{

namespace model
{

namespace affinity
{

constexpr int nwt=1, nwp=1, nl=2, nv=2, nvars=nwt+nwp+nl+nv;
constexpr int ndims=2;
constexpr int nmtx=3;

template<typename Scalar>
struct Affinity
{

  enum
  {
    InputsAtCompileTime=affinity::nvars,
    ValuesAtCompileTime=affinity::ndims
  };

  typedef Scalar T;
  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;

  Affinity(void)
  {}

  static constexpr int
  nvars(void)
  {
    return affinity::nvars;
  }
  static constexpr int
  ndims(void)
  {
    return affinity::ndims;
  }
  static constexpr int
  nmtx(void)
  {
    return affinity::nmtx;
  }

  template<typename U>
  void
  operator()(
    const U* const vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    const Eigen::Map<const vec<U, affinity::nvars> > varsMap(vars);
    const vec<U, affinity::nvars> varst(t*varsMap);
    const Eigen::Map<const vec<U, affinity::nwt> > wtt(varst.template segment<affinity::nwt>(0).data());
    const Eigen::Map<const vec<U, affinity::nwp> > wpt(varst.template segment<affinity::nwp>(affinity::nwt).data());
    const Eigen::Map<const vec<U, affinity::nl> > lt(varst.template segment<nl>(affinity::nwt+affinity::nwp).data());
    const Eigen::Map<const vec<U, affinity::nv> > vt(varst.template segment<nv>(affinity::nwt+affinity::nwp+affinity::nl).data());

    const U costheta=cos(wtt(0)), sintheta=sin(wtt(0));
    const U cosphi=cos(wpt(0)), sinphi=sin(wpt(0));

    mtx<U, 2, 2> rthetaMtx;
    rthetaMtx<<costheta, -sintheta,
      sintheta, costheta;
    mtx<U, 2, 2> rphiMtx;
    rphiMtx<<cosphi, -sinphi,
      sinphi, cosphi;
    mtx<U, 2, 2> dMtx;
    dMtx<<lt.maxCoeff()+U(1.0), U(0.0),
      U(0.0), lt.minCoeff()+U(1.0);

    const mtx<U, 2, 2> aMtx(rthetaMtx*rphiMtx.transpose()*dMtx*rphiMtx);
    m<<aMtx(0, 0), aMtx(0, 1), vt(0),
      aMtx(1, 0), aMtx(1, 1), vt(1),
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

    const Eigen::Map<const vec<V, affinity::ndims> > cMap(c);
    Eigen::Map<vec<U, affinity::ndims> > cmMap(cm);
    vec<U, affinity::ndims+1> ch;
    vec<U, affinity::ndims+1> chm;
    ch<<cMap, V(1.0);
    chm.noalias()=tMtx.inverse()*ch;
    cmMap=chm.template head<affinity::ndims>()/chm(affinity::ndims);
  }

};

}

template<typename T>
using Affinity=Model<affinity::Affinity<T> >;

}

}

#endif // AFFINITY_H
