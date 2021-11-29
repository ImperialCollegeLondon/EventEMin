#ifndef EVENT_EMIN_SIMILARITY_H
#define EVENT_EMIN_SIMILARITY_H

#include "EventEMin/model/model.h"

namespace EventEMin
{
namespace batch
{
template <typename Scalar>
struct Similarity
{
  enum
  {
    NW = 1,
    NS = 1,
    NV = 2,
    NVars = NW + NS + NV,
    NDims = 2,
    NMatrix = 3,
    InputsAtCompileTime = NVars,
    ValuesAtCompileTime = NDims
  };

  typedef Scalar T;
  typedef Vector<T, InputsAtCompileTime> InputType;
  typedef Vector<T, ValuesAtCompileTime> ValueType;

  Similarity(void) = default;

  template <typename U>
  void
  operator()(const U* const vars, const T& t,
             Matrix<U, NMatrix, NMatrix>& m) const
  {
    const Map<const Vector<U, NVars> > varsMap(vars);
    const Vector<U, NVars> varst(t * varsMap);
    const Map<const Vector<U, NW> > wt(varst.template segment<NW>(0).data());
    const Map<const Vector<U, NS> > st(varst.template segment<NS>(NW).data());
    const Map<const Vector<U, NV> > vt(
        varst.template segment<NV>(NW + NS).data());

    const U cosTheta = cos(wt(0)), sinTheta = sin(wt(0));

    Matrix<U, 2, 2> rThetaMatrix;
    rThetaMatrix << cosTheta, -sinTheta, sinTheta, cosTheta;

    const Matrix<U, 2, 2> sMatrix((st(0) + T(1.0)) * rThetaMatrix);
    m << sMatrix(0, 0), sMatrix(0, 1), vt(0), sMatrix(1, 0), sMatrix(1, 1),
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
    Vector<V, NDims + 1> ch;
    Vector<U, NDims + 1> chm;
    ch << cMap, T(1.0);
    chm.noalias() = tMatrix.inverse() * ch;
    cmMap = chm.template head<NDims>() / chm(NDims);
  }
};
}  // namespace batch

template <typename T>
using Similarity = batch::Model<batch::Similarity<T> >;
}  // namespace EventEMin

#endif  // EVENT_EMIN_SIMILARITY_H
