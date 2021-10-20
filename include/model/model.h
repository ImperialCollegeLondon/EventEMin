#ifndef MODEL_H
#define MODEL_H

#include <Eigen/Geometry>
#include <cassert>
#include <cmath>

#include "types_def.h"

namespace event_model
{
template <typename F>
class Model
{
 public:
  typedef F Func;
  typedef typename Func::T T;

  enum
  {
    NVars = Func::NVars,
    NDims = Func::NDims,
    NMatrix = Func::NMatrix,
    InputsAtCompileTime = NVars,
    ValuesAtCompileTime = NDims
  };

  typedef Vector<T, InputsAtCompileTime> InputType;
  typedef Vector<T, ValuesAtCompileTime> ValueType;

 protected:
  const Func func_;

 public:
  Model(void) = default;

  template <typename U>
  void
  operator()(const U* const vars, const T& t,
             Matrix<U, NMatrix, NMatrix>& m) const
  {
    func_(vars, t, m);
  }

  template <typename U>
  void
  operator()(const Vector<U, NVars>& vars, const T& t,
             Matrix<U, NMatrix, NMatrix>& m) const
  {
    (*this)(vars.data(), t, m);
  }

  template <typename U, typename V>
  void
  operator()(const U* const vars, const V* const c, const T& t, U* cm) const
  {
    func_(vars, c, t, cm);
  }

  template <typename U>
  void
  operator()(const Vector<U, NVars>& vars, Vector<U, NDims>* cm,
             const Ref<const Vector<T> >& c, const T& t) const
  {
    assert(c.size() == NDims);

    (*this)(vars.data(), c.data(), t, cm->data());
  }

  template <typename U>
  void
  operator()(const Vector<U, NVars>& vars, const Ref<const Matrix<T> >& c,
             const Ref<const Vector<T> >& ts, const T& tsRef,
             Matrix<U>& cm) const
  {
    assert(c.rows() == NDims);
    assert(cm.rows() == NDims);
    assert(c.cols() == cm.cols());
    assert(c.cols() == ts.size());

    const Vector<T> t(ts.array() - tsRef);
    Vector<U, NDims> vcm;
    for (int k = 0; k < c.cols(); ++k)
    {
      (*this)(vars, &vcm, c.col(k), t(k));
      cm.col(k) = vcm;
    }
  }

  template <typename U>
  void
  operator()(const Vector<U, NVars>& vars, const Ref<const Vector<T> >& c,
             const T& t, Vector<U, NDims>& cm) const
  {
    (*this)(vars, &cm, c, t);
  }

  template <typename U>
  void
  operator()(const Vector<U, NVars>& vars, const Ref<const Matrix<T> >& c,
             const Ref<const Vector<T> >& ts, Matrix<U>& cm) const
  {
    (*this)(vars, c, ts, ts(0), cm);
  }
};
}  // namespace event_model

#endif  // MODEL_H
