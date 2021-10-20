#ifndef TRANSLATION_3D_H
#define TRANSLATION_3D_H

#include "model/model.h"

namespace event_model
{
namespace model
{
template <typename Scalar>
struct Translation3D
{
  enum
  {
    NV = 3,
    NVars = NV,
    NDims = 3,
    NMatrix = 4,
    InputsAtCompileTime = NVars,
    ValuesAtCompileTime = NDims
  };

  typedef Scalar T;
  typedef Vector<T, InputsAtCompileTime> InputType;
  typedef Vector<T, ValuesAtCompileTime> ValueType;

  Translation3D(void) = default;

  template <typename U>
  void
  operator()(const U* const vars, const T& t,
             Matrix<U, NMatrix, NMatrix>& m) const
  {
    const Map<const Vector<U, NVars> > varsMap(vars);
    const Map<const Vector<U, NV> > v(varsMap.template segment<NV>(0).data());

    Matrix<U, NMatrix, NMatrix> skewMtx(Matrix<U, NMatrix, NMatrix>::Zero());
    skewMtx.template rightCols<1>().template head<NV>() = v;
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
}  // namespace model

template <typename T>
using Translation3D = Model<model::Translation3D<T> >;
}  // namespace event_model

#endif  // TRANSLATION_3D_H