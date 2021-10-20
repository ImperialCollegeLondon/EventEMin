#ifndef ROTATION_H
#define ROTATION_H

#include "model/model.h"

namespace event_model
{
namespace model
{
template <typename Scalar>
struct Rotation
{
  enum
  {
    NW = 3,
    NVars = NW,
    NDims = 2,
    NMatrix = 3,
    InputsAtCompileTime = NVars,
    ValuesAtCompileTime = NDims
  };

  typedef Scalar T;
  typedef Vector<T, InputsAtCompileTime> InputType;
  typedef Vector<T, ValuesAtCompileTime> ValueType;

  Rotation(void) = default;

  template <typename U>
  void
  operator()(const U* const vars, const T& t,
             Matrix<U, NMatrix, NMatrix>& m) const
  {
    const Map<const Vector<U, NW> > w(vars);

    m.setIdentity();
    const U wNorm = w.norm();
    if (U(0.0) < wNorm)
    {
      const U theta = t * wNorm;
      Matrix<U, NMatrix, NMatrix> skewMatrix;
      skewMatrix << T(0.0), -w(2), w(1), w(2), T(0.0), -w(0), -w(1), w(0),
          T(0.0);
      skewMatrix /= wNorm;
      m += sin(theta) * skewMatrix +
           (T(1.0) - cos(theta)) * skewMatrix * skewMatrix;
    }
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
    chm.noalias() = tMatrix.transpose() * ch;
    cmMap = chm.template head<NDims>() / chm(NDims);
  }
};
}  // namespace model

template <typename T>
using Rotation = Model<model::Rotation<T> >;
}  // namespace event_model

#endif  // ROTATION_H
