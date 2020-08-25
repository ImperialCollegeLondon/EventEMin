#ifndef SIX_DOF_H
#define SIX_DOF_H

#include "model/model.h"

namespace event_model
{

namespace model
{

namespace six_dof
{

constexpr int nw=3, nv=3, nvars=nw+nv;
constexpr int ndims=3;
constexpr int nmtx=4;

template<typename Scalar>
struct SixDOF
{

  enum
  {
    InputsAtCompileTime=six_dof::nvars,
    ValuesAtCompileTime=six_dof::ndims
  };

  typedef Scalar T;
  typedef vec<T, InputsAtCompileTime> InputType;
  typedef vec<T, ValuesAtCompileTime> ValueType;

  SixDOF(void)
  {}

  static constexpr int
  nvars(void)
  {
    return six_dof::nvars;
  }
  static constexpr int
  ndims(void)
  {
    return six_dof::ndims;
  }
  static constexpr int
  nmtx(void)
  {
    return six_dof::nmtx;
  }

  template<typename U>
  void
  operator()(
    const U* const vars,
    const T& t,
    mtx<U, nmtx(), nmtx()>& m) const
  {
    const Eigen::Map<const vec<U, six_dof::nvars> > varsMap(vars);
    const Eigen::Map<const vec<U, six_dof::nw> > w(varsMap.template segment<six_dof::nw>(0).data());
    const Eigen::Map<const vec<U, six_dof::nw> > v(varsMap.template segment<six_dof::nv>(six_dof::nw).data());

    mtx<U, nmtx(), nmtx()> skewMtx;
    skewMtx<<U(0.0), -w(2), w(1), v(0),
      w(2), U(0.0), -w(0), v(1),
      -w(1), w(0), U(0.0), v(2),
      U(0.0), U(0.0), U(0.0), U(0.0);
    m=mtx<U, nmtx(), nmtx()>::Identity()+t*skewMtx;
    const U wNorm=w.norm();
    if(U(0.0)<wNorm)
    {
      const U theta=t*wNorm;
      const U ctheta=U(1.0)-cos(theta);
      const U stheta=theta-sin(theta);
      skewMtx/=wNorm;
      m+=(ctheta*skewMtx+stheta*skewMtx*skewMtx)*skewMtx;
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

    const Eigen::Map<const vec<V, six_dof::ndims> > cMap(c);
    Eigen::Map<vec<U, six_dof::ndims> > cmMap(cm);
    vec<V, six_dof::ndims+1> ch;
    vec<U, six_dof::ndims+1> chm;
    ch<<cMap, V(1.0);
    chm.noalias()=tMtx.inverse()*ch;
    cmMap=chm.template head<six_dof::ndims>()/chm(six_dof::ndims);
  }

};

}

template<typename T>
using SixDOF=Model<six_dof::SixDOF<T> >;

}

}

#endif // SIX_DOF_H
