#ifndef EVENT_EMIN_HOMOGRAPHY_H
#define EVENT_EMIN_HOMOGRAPHY_H

#include "EventEMin/model/model.h"

namespace EventEMin
{
namespace batch
{
template <typename Scalar>
struct Homography
{
  enum
  {
    NW = 3,
    NV = 3,
    NNSph = 2,
    NN = 3,
    NVars = NW + NV + NNSph,
    NDims = 2,
    NMatrix = 3,
    InputsAtCompileTime = NVars,
    ValuesAtCompileTime = NDims
  };

  typedef Scalar T;
  typedef Vector<T, InputsAtCompileTime> InputType;
  typedef Vector<T, ValuesAtCompileTime> ValueType;

  Homography(void) = default;

  template <typename U>
  void
  operator()(const U* const vars, const T& t,
             Matrix<U, NMatrix, NMatrix>& m) const
  {
    const Map<const Vector<U, NVars> > varsMap(vars);
    const Map<const Vector<U, NW> > w(varsMap.template segment<NW>(0).data());
    const Map<const Vector<U, NV> > v(varsMap.template segment<NV>(NW).data());
    const Map<const Vector<U, NNSph> > nsph(
        varsMap.template segment<NNSph>(NW + NV).data());

    const U cosTheta = cos(nsph(0)), sinTheta = sin(nsph(0));
    const U cosPhi = cos(nsph(1)), sinPhi = sin(nsph(1));
    Vector<U, NN> n;
    n << cosTheta * sinPhi, sinTheta * sinPhi, cosPhi;

    Matrix<U, NMatrix, NMatrix> rMatrix(
        Matrix<U, NMatrix, NMatrix>::Identity());
    const U wNorm = w.norm();
    if (U(0.0) < wNorm)
    {
      const U theta = t * wNorm;
      Matrix<U, NMatrix, NMatrix> skewMatrix;
      skewMatrix << T(0.0), -w(2), w(1), w(2), T(0.0), -w(0), -w(1), w(0),
          T(0.0);
      skewMatrix /= wNorm;
      rMatrix += sin(theta) * skewMatrix +
                 (T(1.0) - cos(theta)) * skewMatrix * skewMatrix;
    }
    m = rMatrix - t * v * n.transpose();
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
using Homography = batch::Model<batch::Homography<T> >;
}  // namespace EventEMin

#endif  // EVENT_EMIN_HOMOGRAPHY_H
