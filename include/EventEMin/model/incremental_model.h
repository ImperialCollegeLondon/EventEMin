#ifndef EVENT_EMIN_INCREMENTAL_MODEL_H
#define EVENT_EMIN_INCREMENTAL_MODEL_H

#include <cassert>
#include <cmath>

#include "EventEMin/types_def.h"

namespace EventEMin
{
namespace incremental
{
template <typename F>
class Model
{
 public:
  typedef F Func;
  typedef typename Func::T T;
  typedef typename Func::Point Point;
  typedef typename Func::PointHomogeneous PointHomogeneous;
  typedef typename Func::Vars Vars;
  typedef typename Func::CMatrix CMatrix;
  typedef typename Func::DMatrix DMatrix;
  typedef typename Func::PMatrix PMatrix;
  typedef typename Func::GMatrix GMatrix;
  typedef typename Func::TMatrix TMatrix;

  enum
  {
    NVars = Func::NVars,
    NDims = Func::NDims,
    NMatrix = Func::NMatrix,
  };

 protected:
  Func func_;

 public:
  Model(void) = default;

  void
  operator()(const Vars& vars, const Ref<const Matrix<T> >& c,
             const Ref<const Vector<T> >& ts, Matrix<T>& cm) const
  {
    assert(c.rows() == NDims);
    assert(cm.rows() == NDims);
    assert(c.cols() == cm.cols());
    assert(c.cols() == ts.size());

    const Vector<T> t(ts.array() - ts(0));
    for (int k = 0; k < c.cols(); ++k)
    {
      (*this)(vars, c.col(k), t(k), cm.col(k));
    }
  }

  void
  operator()(const Ref<const Vars>& vars, const Ref<const Point>& c, const T& t,
             Ref<Point> cm) const
  {
    func_(vars, c, t, cm);
  }

  void
  operator()(const Ref<const Vars>& vars, const Ref<const Point>& c, const T& t,
             Ref<Point> cm, Ref<DMatrix> dcm, Ref<Point> cg,
             Ref<PMatrix> pMtx) const
  {
    func_(vars, c, t, cm, dcm, cg, pMtx);
  }

  void
  operator()(const Ref<const Vars>& vars, const Ref<const Point>& c, const T& t,
             Ref<Point> cm, Ref<DMatrix> dcm, Ref<CMatrix> dc) const
  {
    func_(vars, c, t, cm, dcm, dc);
  }

  void
  transformation(const Ref<const Vars>& varst, Ref<TMatrix> tMatrix) const
  {
    func_.transformation(varst, tMatrix);
  }
};
}  // namespace incremental
}  // namespace EventEMin

#endif  // EVENT_EMIN_INCREMENTAL_MODEL_H
