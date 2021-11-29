#ifndef EVENT_EMIN_INCREMENTAL_SIX_DOF_H
#define EVENT_EMIN_INCREMENTAL_SIX_DOF_H

#include "EventEMin/model/incremental_model.h"

namespace EventEMin
{
namespace incremental
{
template <typename Scalar>
struct SixDOF
{
  enum
  {
    NW = 3,
    NV = 3,
    NVars = NW + NV,
    NDims = 3,
    NMatrix = 4
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

  SixDOF(void) = default;

  void
  operator()(const Ref<const Vars>& vars, const Ref<const Point>& c, const T& t,
             Ref<Point> cm) const
  {
    const Vars varst(t * vars);

    GMatrix gMatrix;
    generator(varst, gMatrix);

    PointHomogeneous ch;
    ch << c, T(1.0);

    // transformation
    transformation(gMatrix, ch, cm);
  }

  void
  operator()(const Ref<const Vars>& vars, const Ref<const Point>& c, const T& t,
             Ref<Point> cm, Ref<DMatrix> dcm, Ref<Point> cg,
             Ref<PMatrix> pMatrix) const
  {
    const Vars varst(t * vars);

    GMatrix gMatrix;
    generator(varst, gMatrix);

    PointHomogeneous ch;
    ch << c, T(1.0);

    // transformation
    transformation(gMatrix, ch, cm);

    // generator
    generator(gMatrix, ch, cg);

    // derivatives wrt to vars
    dtransformationdvars(ch, t, dcm);

    // perturbation
    perturbation(c, t, pMatrix);
  }

  void
  operator()(const Ref<const Vars>& vars, const Ref<const Point>& c, const T& t,
             Ref<Point> cm, Ref<DMatrix> dcm, Ref<CMatrix> dc) const
  {
    const Vars varst(t * vars);

    GMatrix gMatrix;
    generator(varst, gMatrix);

    PointHomogeneous ch;
    ch << c, T(1.0);

    // transformation
    TMatrix tMatrix;
    transformation(gMatrix, ch, cm, tMatrix);

    // derivatives wrt to vars
    dtransformationdvars(ch, t, dcm);

    // derivatives wrt to c
    dc = tMatrix.template topLeftCorner<NDims, NDims>();
  }

  void
  generator(const Ref<const GMatrix>& gMatrix,
            const Ref<const PointHomogeneous>& ch, Ref<Point> cg) const
  {
    cg.noalias() = gMatrix.template topRows<NDims>() * ch;
  }

  void
  generator(const Ref<const Vars>& varst, Ref<GMatrix> gMatrix) const
  {
    gMatrix << T(0.0), -varst(2), varst(1), varst(3), varst(2), T(0.0),
        -varst(0), varst(4), -varst(1), varst(0), T(0.0), varst(5), T(0.0),
        T(0.0), T(0.0), T(0.0);
  }

  void
  dgeneratordwx(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), -t, T(0.0),
        T(0.0), t, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0);
  }

  void
  dgeneratordwy(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), T(0.0), t, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), -t,
        T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0);
  }

  void
  dgeneratordwz(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), -t, T(0.0), T(0.0), t, T(0.0), T(0.0), T(0.0), T(0.0),
        T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0);
  }

  void
  dgeneratordvx(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), T(0.0), T(0.0), t, T(0.0), T(0.0), T(0.0), T(0.0),
        T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0);
  }

  void
  dgeneratordvy(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), t,
        T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0);
  }

  void
  dgeneratordvz(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
        T(0.0), T(0.0), T(0.0), t, T(0.0), T(0.0), T(0.0), T(0.0);
  }

  void
  perturbation(const Ref<const Point>& c, const T& t,
               Ref<PMatrix> pMatrix) const
  {
    const Point ct(t * c);
    pMatrix << T(0.0), ct(2), -ct(1), t, T(0.0), T(0.0), -ct(2), T(0.0), ct(0),
        T(0.0), t, T(0.0), ct(1), -ct(0), T(0.0), T(0.0), T(0.0), t;
  }

  void
  transformation(const Ref<const GMatrix>& gMatrix,
                 const Ref<const PointHomogeneous>& ch, Ref<Point> cm) const
  {
    TMatrix tMatrix;
    transformation(gMatrix, ch, cm, tMatrix);
  }

  void
  transformation(const Ref<const GMatrix>& gMatrix,
                 const Ref<const PointHomogeneous>& ch, Ref<Point> cm,
                 Ref<TMatrix> tMatrix) const
  {
    transformationG(gMatrix, tMatrix);
    cm.noalias() = tMatrix.template topRows<NDims>() * ch;
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
    // 1st order Taylor series expansion
    tMatrix.setIdentity();
    tMatrix += gMatrix;
  }

  void
  dtransformationdvars(const Ref<const PointHomogeneous>& ch, const T& t,
                       Ref<DMatrix> dcm) const
  {
    GMatrix dgMatrix;
    dgeneratordwx(t, dgMatrix);
    dtransformationdvars(dgMatrix, ch, dcm.col(0));
    dgeneratordwy(t, dgMatrix);
    dtransformationdvars(dgMatrix, ch, dcm.col(1));
    dgeneratordwz(t, dgMatrix);
    dtransformationdvars(dgMatrix, ch, dcm.col(2));
    dgeneratordvx(t, dgMatrix);
    dtransformationdvars(dgMatrix, ch, dcm.col(3));
    dgeneratordvy(t, dgMatrix);
    dtransformationdvars(dgMatrix, ch, dcm.col(4));
    dgeneratordvz(t, dgMatrix);
    dtransformationdvars(dgMatrix, ch, dcm.col(5));
  }

  void
  dtransformationdvars(const Ref<const GMatrix>& dgMatrix,
                       const Ref<const PointHomogeneous>& ch,
                       Ref<Point> dcm) const
  {
    dcm.noalias() = dgMatrix.template topRows<NDims>() * ch;
  }
};
}  // namespace incremental

template <typename T>
using IncrementalSixDOF = incremental::Model<incremental::SixDOF<T> >;
}  // namespace EventEMin

#endif  // EVENT_EMIN_INCREMENTAL_SIX_DOF_H
