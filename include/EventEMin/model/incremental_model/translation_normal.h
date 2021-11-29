#ifndef EVENT_EMIN_INCREMENTAL_TRANSLATION_NORMAL_H
#define EVENT_EMIN_INCREMENTAL_TRANSLATION_NORMAL_H

#include "EventEMin/model/incremental_model.h"

namespace EventEMin
{
namespace incremental
{
template <typename Scalar>
struct TranslationNormal
{
  enum
  {
    NV = 3,
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

  TranslationNormal(void) = default;

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
    const T z = transformation(varst, c, cm);

    PointHomogeneous ch;
    ch << c, T(1.0);

    // generator
    generator(varst, ch, z, cg);

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

    // transformation
    const T z = transformation(varst, c, cm);

    PointHomogeneous ch;
    ch << c, T(1.0);

    // derivatives wrt to vars
    dtransformationdvars(ch, cm, t, z, dcm);

    // derivatives wrt to c
    const T zInv = T(1.0) / z;
    dc << zInv, T(0.0), T(0.0), zInv;
  }

  void
  generator(const Ref<const Vars>& varst, const Ref<const PointHomogeneous>& ch,
            const T& z, Ref<Point> cg) const
  {
    GMatrix gMatrix;
    generator(varst, gMatrix);
    cg.noalias() = (gMatrix.template topRows<NDims>() * ch) / z;
  }

  void
  generator(const Ref<const Vars>& varst, Ref<GMatrix> gMatrix) const
  {
    gMatrix << varst(2), T(0.0), -varst(0), T(0.0), varst(2), -varst(1), T(0.0),
        T(0.0), T(0.0);
  }

  void
  perturbation(const Ref<const Point>& c, const T& t, const T& z,
               Ref<PMatrix> pMatrix) const
  {
    const T tz = t / z;
    const Point ctz(tz * c);
    pMatrix << -tz, T(0.0), ctz(0), T(0.0), -tz, ctz(1);
  }

  T
  transformation(const Ref<const Vars>& varst, const Ref<const Point>& c,
                 Ref<Point> cm) const
  {
    const T z = T(1.0) - varst(2);
    cm.noalias() = (c - varst.template head<2>()) / z;
    return z;
  }

  void
  dtransformationdvx(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), T(0.0), -t, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
        T(0.0);
  }

  void
  dtransformationdvy(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), -t, T(0.0), T(0.0),
        T(0.0);
  }

  void
  dtransformationdvz(const T& t, Ref<GMatrix> dgMatrix) const
  {
    dgMatrix << T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
        -t;
  }

  void
  dtransformationdvars(const Ref<const PointHomogeneous>& ch,
                       const Ref<const Point>& cm, const T& t, const T& z,
                       Ref<DMatrix> dcm) const
  {
    GMatrix dgMatrix;
    dtransformationdvx(t, dgMatrix);
    dtransformationdvars(dgMatrix, ch, cm, z, dcm.col(0));
    dtransformationdvy(t, dgMatrix);
    dtransformationdvars(dgMatrix, ch, cm, z, dcm.col(1));
    dtransformationdvz(t, dgMatrix);
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
};
}  // namespace incremental

template <typename T>
using IncrementalTranslationNormal =
    incremental::Model<incremental::TranslationNormal<T> >;
}  // namespace EventEMin

#endif  // EVENT_EMIN_INCREMENTAL_TRANSLATION_NORMAL_H
