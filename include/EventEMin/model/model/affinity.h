#ifndef EVENT_EMIN_AFFINITY_H
#define EVENT_EMIN_AFFINITY_H

#include "EventEMin/model/model.h"

namespace EventEMin
{
namespace batch
{
template <typename Scalar>
struct Affinity
{
  enum
  {
    NWt = 1,
    NWp = 1,
    NL = 2,
    NV = 2,
    NVars = NWt + NWp + NL + NV,
    NDims = 2,
    NMatrix = 3,
    InputsAtCompileTime = NVars,
    ValuesAtCompileTime = NDims
  };

  typedef Scalar T;
  typedef Vector<T, InputsAtCompileTime> InputType;
  typedef Vector<T, ValuesAtCompileTime> ValueType;

  Affinity(void) = default;

  template <typename U>
  void
  operator()(const U* const vars, const T& t,
             Matrix<U, NMatrix, NMatrix>& m) const
  {
    const Map<const Vector<U, NVars> > varsMap(vars);
    const Vector<U, NVars> varst(t * varsMap);
    const Map<const Vector<U, NWt> > wtt(varst.template segment<NWt>(0).data());
    const Map<const Vector<U, NWp> > wpt(
        varst.template segment<NWp>(NWt).data());
    const Map<const Vector<U, NL> > lt(
        varst.template segment<NL>(NWt + NWp).data());
    const Map<const Vector<U, NV> > vt(
        varst.template segment<NV>(NWt + NWp + NL).data());

    const U cosTheta = cos(wtt(0)), sinTheta = sin(wtt(0));
    const U cosPhi = cos(wpt(0)), sinPhi = sin(wpt(0));

    Matrix<U, 2, 2> rThetaMatrix;
    rThetaMatrix << cosTheta, -sinTheta, sinTheta, cosTheta;
    Matrix<U, 2, 2> rPhiMatrix;
    rPhiMatrix << cosPhi, -sinPhi, sinPhi, cosPhi;
    Matrix<U, 2, 2> dMatrix;
    dMatrix << lt.maxCoeff() + T(1.0), T(0.0), T(0.0), lt.minCoeff() + T(1.0);

    const Matrix<U, 2, 2> aMatrix(rThetaMatrix * rPhiMatrix.transpose() *
                                  dMatrix * rPhiMatrix);
    m << aMatrix(0, 0), aMatrix(0, 1), vt(0), aMatrix(1, 0), aMatrix(1, 1),
        vt(1), T(0.0), T(0.0), T(1.0);
  }

  template <typename U, typename V>
  void
  operator()(const U* const vars, const V* const c, const T& t, U* cm) const
  {
    Matrix<U, NMatrix, NMatrix> tMatrix;
    (*this)(vars, t, tMatrix);

    const Map<const Vector<V, NDims> > cMap(c);
    Map<Vector<U, NDims> > cmMap(cm);
    Vector<U, NDims + 1> ch;
    Vector<U, NDims + 1> chm;
    ch << cMap, T(1.0);
    chm.noalias() = tMatrix.inverse() * ch;
    cmMap = chm.template head<NDims>() / chm(NDims);
  }
};
}  // namespace batch

template <typename T>
using Affinity = batch::Model<batch::Affinity<T> >;
}  // namespace EventEMin

#endif  // EVENT_EMIN_AFFINITY_H
