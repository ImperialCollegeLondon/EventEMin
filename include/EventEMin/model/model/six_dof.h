#ifndef EVENT_EMIN_SIX_DOF_H
#define EVENT_EMIN_SIX_DOF_H

#include "EventEMin/model/model.h"

namespace EventEMin
{
namespace batch
{
template <typename Scalar>
struct SixDOF
{
  enum
  {
    NW = 3,
    NV = 3,
    NVars = NW + NV,
    NDims = 3,
    NMatrix = 4,
    InputsAtCompileTime = NVars,
    ValuesAtCompileTime = NDims
  };

  typedef Scalar T;
  typedef Vector<T, InputsAtCompileTime> InputType;
  typedef Vector<T, ValuesAtCompileTime> ValueType;

  SixDOF(void) = default;

  template <typename U>
  void
  operator()(const U* const vars, const T& t,
             Matrix<U, NMatrix, NMatrix>& m) const
  {
    const Map<const Vector<U, NVars> > varsMap(vars);
    const Map<const Vector<U, NW> > w(varsMap.template segment<NW>(0).data());
    const Map<const Vector<U, NV> > v(varsMap.template segment<NV>(NW).data());

    // First-order Taylor expansion
    Matrix<U, NMatrix, NMatrix> skewMtx;
    skewMtx << T(0.0), -w(2), w(1), v(0), w(2), T(0.0), -w(0), v(1), -w(1),
        w(0), T(0.0), v(2), T(0.0), T(0.0), T(0.0), T(0.0);
    m = t * skewMtx;
    m.diagonal().array() += T(1.0);
  }

  template <typename U, typename V>
  void
  operator()(const U* const vars, const V* const c, const T& t, U* cm) const
  {
    Matrix<U, NMatrix, NMatrix> tMtx;
    (*this)(vars, t, tMtx);

    const Map<const Vector<V, NDims> > cMap(c);
    Map<Vector<U, NDims> > cmMap(cm);
    Vector<V, NDims + 1> ch;
    Vector<U, NDims + 1> chm;
    ch << cMap, T(1.0);
    chm.noalias() = tMtx.inverse() * ch;
    cmMap = chm.template head<NDims>() / chm(NDims);
  }
};
}  // namespace batch

template <typename T>
using SixDOF = batch::Model<batch::SixDOF<T> >;
}  // namespace EventEMin

#endif  // EVENT_EMIN_SIX_DOF_H
