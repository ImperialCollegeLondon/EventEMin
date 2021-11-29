#ifndef EVENT_EMIN_INCREMENTAL_ROTATION_H
#define EVENT_EMIN_INCREMENTAL_ROTATION_H

#include "EventEMin/model/incremental_model.h"

namespace EventEMin
{
namespace incremental
{
template <typename Scalar>
struct Rotation
{
  enum
  {
    NW = 3,
    NVars = NW,
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

  Rotation(void) = default;

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
    const T z = transformation(gMatrix, ch, cm);

    // generator
    generator(gMatrix, ch, z, cg);

    // derivatives wrt to vars
    dtransformationdvars(ch, cm, t, z, dcm);

    // perturbation
    perturbation(c, t, z, pMatrix);
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
    const T z = transformation(gMatrix, ch, cm, tMatrix);

    // derivatives wrt to vars
    dtransformationdvars(ch, cm, t, z, dcm);

    // derivatives wrt to c
    dtransformationdc(tMatrix, z, dc);
  }

  void
  generator(const Ref<const GMatrix>& gMatrix,
            const Ref<const PointHomogeneous>& ch, const T& z,
            Ref<Point> cg) const
  {
    cg.noalias() = (gMatrix.template topRows<NDims>() * ch) / z;
  }

  void
  generator(const Ref<const Vars>& varst, Ref<GMatrix> gMatrix) const
  {
    gMatrix << T(0.0), -varst(2), varst(1), varst(2), T(0.0), -varst(0),
        -varst(1), varst(0), T(0.0);
  }

  void
  dgeneratordwx(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), -t, T(0.0), t, T(0.0);
  }

  void
  dgeneratordwy(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), T(0.0), t, T(0.0), T(0.0), T(0.0), -t, T(0.0), T(0.0);
  }

  void
  dgeneratordwz(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), -t, T(0.0), t, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0);
  }

  void
  perturbation(const Ref<const Point>& c, const T& t, const T& z,
               Ref<PMatrix> pMatrix) const
  {
    const T tz = t / z;
    const Point ctz(tz * c);
    pMatrix << T(0.0), tz, -ctz(1), -tz, T(0.0), ctz(0);
  }

  T
  transformation(const Ref<const GMatrix>& gMatrix,
                 const Ref<const PointHomogeneous>& ch, Ref<Point> cm) const
  {
    TMatrix tMatrix;
    return transformation(gMatrix, ch, cm, tMatrix);
  }

  T
  transformation(const Ref<const GMatrix>& gMatrix,
                 const Ref<const PointHomogeneous>& ch, Ref<Point> cm,
                 Ref<TMatrix> tMatrix) const
  {
    transformationG(gMatrix, tMatrix);
    const PointHomogeneous cmh(tMatrix * ch);
    cm = cmh.template head<NDims>() / cmh(2);
    return cmh(2);
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
  dtransformationdvars(const Ref<const PointHomogeneous>& ch,
                       const Ref<const Point>& cm, const T& t, const T& z,
                       Ref<DMatrix> dcm) const
  {
    GMatrix dgMatrix;
    dgeneratordwx(t, dgMatrix);
    dtransformationdvars(dgMatrix, ch, cm, z, dcm.col(0));
    dgeneratordwy(t, dgMatrix);
    dtransformationdvars(dgMatrix, ch, cm, z, dcm.col(1));
    dgeneratordwz(t, dgMatrix);
    dtransformationdvars(dgMatrix, ch, cm, z, dcm.col(2));
  }

  void
  dtransformationdvars(const Ref<const GMatrix>& dgMatrix,
                       const Ref<const PointHomogeneous>& ch,
                       const Ref<const Point>& cm, const T& z,
                       Ref<Point> dcm) const
  {
    dcm = (dgMatrix.template topRows<NDims>() * ch -
           dgMatrix.template bottomRows<1>() * ch * cm) /
          z;
  }

  void
  dtransformationdc(const Ref<const TMatrix>& tMatrix, const T& z,
                    Ref<CMatrix> dc) const
  {
    PointHomogeneous dcc;
    dcc << T(1.0), T(0.0), T(0.0);
    dtransformationdc(tMatrix, dcc, z, dc.col(0));
    dcc << T(0.0), T(1.0), T(0.0);
    dtransformationdc(tMatrix, dcc, z, dc.col(1));
  }

  void
  dtransformationdc(const Ref<const TMatrix>& tMatrix,
                    const Ref<const PointHomogeneous>& dcc, const T& z,
                    Ref<Point> dc) const
  {
    dc.noalias() =
        tMatrix.template topRows<NDims>() *
        (dcc.array() - (tMatrix.template bottomRows<1>() * dcc)(0) / z)
            .matrix() /
        z;
  }
};
}  // namespace incremental

template <typename T>
using IncrementalRotation = incremental::Model<incremental::Rotation<T> >;
}  // namespace EventEMin

#endif  // EVENT_EMIN_INCREMENTAL_ROTATION_H
