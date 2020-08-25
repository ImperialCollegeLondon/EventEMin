#ifndef HOMOGRAPHY_H
#define HOMOGRAPHY_H

#include "model/model.h"

namespace event_model
{

namespace model
{

namespace homography
{

constexpr int nw=3, nv=3, nnsph=2, nn=3, nvars=nw+nv+nnsph;
constexpr int ndims=2;
constexpr int nmtx=3;

template<typename Scalar>
struct Homography
{

  enum
  {
    InputsAtCompileTime=homography::nvars,
    ValuesAtCompileTime=homography::ndims
  };

  typedef Scalar T;
  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;

  Homography(void)
  {}

  static constexpr int
  nvars(void)
  {
    return homography::nvars;
  }
  static constexpr int
  ndims(void)
  {
    return homography::ndims;
  }
  static constexpr int
  nmtx(void)
  {
    return homography::nmtx;
  }

  template<typename U>
  void
  operator()(
    const U* const vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    const Eigen::Map<const vec<U, homography::nvars> > varsMap(vars);
    const Eigen::Map<const vec<U, homography::nw> > w(varsMap.template segment<homography::nw>(0).data());
    const Eigen::Map<const vec<U, homography::nv> > v(varsMap.template segment<homography::nv>(homography::nw).data());
    const Eigen::Map<const vec<U, homography::nnsph> > nsph(varsMap.template segment<homography::nnsph>(homography::nw+homography::nv).data());

    const U costheta=cos(nsph(0)), sintheta=sin(nsph(0));
    const U cosphi=cos(nsph(1)), sinphi=sin(nsph(1));
    vec<U, homography::nn> n;
    n<<costheta*sinphi, sintheta*sinphi, cosphi;

    mtx<U, nmtx(), nmtx()> rMtx(mtx<U, nmtx(), nmtx()>::Identity());
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
      rMtx+=stheta*skewMtx+ctheta*skewMtx*skewMtx;
    }
    m=rMtx-t*v*n.transpose();
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

    const Eigen::Map<const vec<V, homography::ndims> > cMap(c);
    Eigen::Map<vec<U, homography::ndims> > cmMap(cm);
    vec<V, homography::ndims+1> ch;
    vec<U, homography::ndims+1> chm;
    ch<<cMap, V(1.0);
    chm.noalias()=tMtx.inverse()*ch;
    cmMap=chm.template head<homography::ndims>()/chm(homography::ndims);
  }

};

}

template<typename T>
using Homography=Model<homography::Homography<T> >;

}

}

#endif // HOMOGRAPHY_H
