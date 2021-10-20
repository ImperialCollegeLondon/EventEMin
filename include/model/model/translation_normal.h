#ifndef TRANSLATION_NORMAL_H
#define TRANSLATION_NORMAL_H

#include "model/model.h"

namespace event_model
{
namespace model
{
template <typename Scalar>
struct TranslationNormal
{
  enum
  {
    NV = 3,
    NVars = NV,
    NDims = 2,
    NMatrix = 3,
    InputsAtCompileTime = NVars,
    ValuesAtCompileTime = NDims
  };

  typedef Scalar T;
  typedef Vector<T, InputsAtCompileTime> InputType;
  typedef Vector<T, ValuesAtCompileTime> ValueType;

  TranslationNormal(void) = default;

  template <typename U>
  void
  operator()(const U* const vars, const T& t,
             Matrix<U, NMatrix, NMatrix>& m) const
  {
    const Map<const Vector<U, NV> > v(vars);
    const Vector<U, NV> vt(t * v);

    m.setIdentity();
    m.template rightCols<1>() -= vt;
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
using TranslationNormal = Model<model::TranslationNormal<T> >;
}  // namespace event_model

#endif  // TRANSLATION_NORMAL_H
