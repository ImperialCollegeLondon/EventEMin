#ifndef EVENT_EMIN_INCREMENTAL_TRANSLATION_2D_H
#define EVENT_EMIN_INCREMENTAL_TRANSLATION_2D_H

#include "EventEMin/model/incremental_model.h"

namespace EventEMin
{
namespace incremental
{
template <typename Scalar>
struct Translation2D
{
  enum
  {
    NV = 2,
    NVars = NV,
    NDims = 2,
    NMatrix = 3
  };

  typedef Scalar T;
  typedef Vector<T, NDims> Point;
  typedef Vector<T, NDims + 1> PointHomogeneous;
  typedef Vector<T, NVars> Vars;
  typedef Matrix<T, NDims, NDims> CMatrix;
  typedef Matrix<T, NDims, NVars> DMatrix;
  typedef Matrix<T, NDims, NVars> PMatrix;
  typedef Matrix<T, NMatrix, NMatrix> GMatrix;
  typedef Matrix<T, NMatrix, NMatrix> TMatrix;

  Translation2D(void) = default;

  void
  operator()(const Ref<const Vars>& vars, const Ref<const Point>& c, const T& t,
             Ref<Point> cm) const
  {
    const Vars varst(t * vars);

    // transformation
    transformation(varst, c, cm);
  }

  void
  operator()(const Ref<const Vars>& vars, const Ref<const Point>& c, const T& t,
             Ref<Point> cm, Ref<DMatrix> dcm, Ref<Point> cg,
             Ref<PMatrix> pMatrix) const
  {
    const Vars varst(t * vars);

    // transformation
    transformation(varst, c, cm);

    // generator
    cg = varst;

    // derivatives wrt to vars
    dcm << t, T(0.0), T(0.0), t;

    // perturbation
    pMatrix << t, T(0.0), T(0.0), t;
  }

  void
  operator()(const Ref<const Vars>& vars, const Ref<const Point>& c, const T& t,
             Ref<Point> cm, Ref<DMatrix> dcm, Ref<CMatrix> dc) const
  {
    const Vars varst(t * vars);

    // transformation
    transformation(varst, c, cm);

    // derivatives wrt to vars
    dcm << t, T(0.0), T(0.0), t;

    // derivatives wrt to c
    dc.setIdentity();
  }

  void
  generator(const Ref<const Vars>& varst, Ref<GMatrix> gMatrix) const
  {
    gMatrix << T(0.0), T(0.0), varst(0), T(0.0), T(0.0), varst(1), T(0.0),
        T(0.0), T(0.0);
  }

  void
  transformation(const Ref<const Vars>& varst, const Ref<const Point>& c,
                 Ref<Point> cm) const
  {
    cm = c + varst;
  }

  void
  transformation(const Ref<const Vars>& varst, Ref<TMatrix> tMatrix) const
  {
    GMatrix gMatrix;
    generator(varst, gMatrix);
    transformationG(gMatrix, tMatrix);
  }

  void
  transformationG(const Ref<const GMatrix>& gMatrix, Ref<TMatrix> tMatrix) const
  {
    tMatrix.setIdentity();
    tMatrix += gMatrix;
  }
};
}  // namespace incremental

template <typename T>
using IncrementalTranslation2D =
    incremental::Model<incremental::Translation2D<T> >;
}  // namespace EventEMin

#endif  // EVENT_EMIN_INCREMENTAL_TRANSLATION_2D_H
