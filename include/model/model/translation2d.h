#ifndef TRANSLATION_2D_H
#define TRANSLATION_2D_H

#include "model/model.h"

namespace event_model
{
namespace model
{
template <typename Scalar>
struct Translation2D
{
  enum
  {
    NV = 2,
    NVars = NV,
    NDims = 2,
    NMatrix = 3,
    InputsAtCompileTime = NVars,
    ValuesAtCompileTime = NDims
  };

  typedef Scalar T;
  typedef Vector<T, InputsAtCompileTime> InputType;
  typedef Vector<T, ValuesAtCompileTime> ValueType;

  Translation2D(void) = default;

  template <typename U>
  void
  operator()(const U* const vars, const T& t,
             Matrix<U, NMatrix, NMatrix>& m) const
  {
    const Map<const Vector<U, NV> > v(vars);

    m.setIdentity();
    m.template topRightCorner<NV, 1>() = t * v;
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
}  // namespace model

template <typename T>
using Translation2D = Model<model::Translation2D<T> >;
}  // namespace event_model

#endif  // TRANSLATION_2D_H
