#ifndef EVENT_EMIN_TRANSLATION_H
#define EVENT_EMIN_TRANSLATION_H

#include "EventEMin/model/model.h"

namespace EventEMin
{
namespace batch
{
template <typename Scalar>
struct Translation
{
  enum
  {
    NV = 3,
    NNSph = 2,
    NN = 3,
    NVars = NV + NNSph,
    NDims = 2,
    NMatrix = 3,
    InputsAtCompileTime = NVars,
    ValuesAtCompileTime = NDims
  };

  typedef Scalar T;
  typedef Vector<T, InputsAtCompileTime> InputType;
  typedef Vector<T, ValuesAtCompileTime> ValueType;

  Translation(void) = default;

  template <typename U>
  void
  operator()(const U* const vars, const T& t,
             Matrix<U, NMatrix, NMatrix>& m) const
  {
    const Map<const Vector<U, NVars> > varsMap(vars);
    const Vector<U, NV> varst(t * varsMap.template head<NV>());
    const Map<const Vector<U, NV> > vt(varst.template segment<NV>(0).data());
    const Map<const Vector<U, NNSph> > nsph(
        varsMap.template segment<NNSph>(NV).data());

    const U cosTheta = cos(nsph(0)), sinTheta = sin(nsph(0));
    const U cosphi = cos(nsph(1)), sinphi = sin(nsph(1));
    Vector<U, NN> n;
    n << cosTheta * sinphi, sinTheta * sinphi, cosphi;

    m.noalias() = -vt * n.transpose();
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
using Translation = batch::Model<batch::Translation<T> >;
}  // namespace EventEMin

#endif  // EVENT_EMIN_TRANSLATION_H
