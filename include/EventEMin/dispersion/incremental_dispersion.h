#ifndef EVENT_EMIN_INCREMENTAL_DISPERSION_H
#define EVENT_EMIN_INCREMENTAL_DISPERSION_H

#include <cassert>

#include "EventEMin/data_stats.h"
#include "EventEMin/types_def.h"
#include "EventEMin/utilities.h"

namespace EventEMin
{
namespace incremental
{
template <typename T>
struct DispersionParams
{
  T minStep;
  int maxIter, wSize;

  DispersionParams(const T& minStep = T(1.0e-6), const int maxIter = 10,
                   const int wSize = 4)
      : minStep(minStep), maxIter(maxIter), wSize(wSize)
  {
  }
};

class Cell
{
 private:
  int n_, nMax_, nMaxMinus1_;
  int first_, last_;

  StdVector<int> ind_;

 public:
  Cell(const int nMax = 2)
      : n_(0),
        nMax_(nMax),
        nMaxMinus1_(nMax_ - 1),
        first_(0),
        last_(0),
        ind_(nMax_)
  {
    assert(0 < nMaxMinus1_ && (nMax_ & nMaxMinus1_) == 0);
  }

  int
  n(void) const
  {
    return n_;
  }
  int
  ind(const int k) const
  {
    return ind_[fastAddCyclic(first_, k, nMaxMinus1_)];
  }

  void
  add(const int ind)
  {
    if (n_ >= nMax_)
    {
      nMax_ <<= 1;
      nMaxMinus1_ = nMax_ - 1;
      ind_.resize(nMax_);
      if (0 < first_)
      {
        std::copy_n(ind_.begin(), last_, ind_.begin() + n_);
      }
      last_ += n_;
    }

    ind_[last_] = ind;
    ++n_;
    last_ = fastIncrementCyclic(last_, nMaxMinus1_);
  }

  int
  remove(void)
  {
    if (0 < n_)
    {
      first_ = fastIncrementCyclic(first_, nMaxMinus1_);
      --n_;
    }
    return n_;
  }
};

template <typename M>
struct Stats
{
  typedef typename M::Point Point;
  typedef typename M::CMatrix CMatrix;

  Point mean, singularValues;
  CMatrix cov, w;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Stats(void)
      : mean(Point::Zero()),
        singularValues(Point::Ones()),
        cov(CMatrix::Identity()),
        w(CMatrix::Identity())
  {
  }

  void
  update(const Ref<const Point>& c, const int n)
  {
    const Point cDiff(c - mean);
    mean += cDiff / n;
    cov.noalias() += (cDiff * (c - mean).transpose() - cov) / n;

    computeWhitening(cov, w, singularValues);
  }
};

template <typename M, bool W = true>
class Whitening
{
 public:
  typedef M Model;
  typedef typename Model::T T;

  typedef typename Model::Point Point;
  typedef typename Model::Vars Vars;

 private:
  projectEvent<T, Model::NDims> projectEvent_;

  Model model_;

  Stats<M> istats_;

 public:
  Whitening(void) = default;

  void
  processEvent(const Matrix<T, 3, 3>& camParams, const Ref<const Point>& scale,
               const Ref<const Point>& c, Ref<Point> cs) const
  {
    if constexpr (W)
    {
      Point cw;
      whitenPoints(c - istats_.mean, istats_.w, cw);
      projectEvent_(
          camParams,
          ((istats_.singularValues.array() * cw.array() + istats_.mean.array())
               .matrix()),
          cs);
    }
    else
    {
      projectEvent_(camParams, c, cs);
    }
    cs.array() *= scale.array();
  }

  void
  updateStats([[maybe_unused]] const Ref<const Vars>& vars,
              [[maybe_unused]] const Ref<const Point>& c,
              [[maybe_unused]] const T& tsDiffRef, [[maybe_unused]] const int n)
  {
    if constexpr (W)
    {
      Point cm;
      model_(vars, c, tsDiffRef, cm);
      istats_.update(cm, n);
    }
  }
};

template <typename M>
struct VarsEstimate
{
  typedef M Model;
  typedef typename Model::T T;

  enum
  {
    NVars = Model::NVars
  };

  typedef typename Model::Vars Vars;
  typedef Matrix<T, NVars, NVars> VarsMatrix;

  T val, valCum;
  Vars vNum, vNumCum;
  VarsMatrix vDen, vDenCum;

  Vars vars;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  VarsEstimate(void)
      : valCum(T(0.0)),
        vNumCum(Vars::Zero()),
        vDenCum(VarsMatrix::Zero()),
        vars(Vars::Zero())
  {
  }

  void
  decay(const T& dec)
  {
    valCum *= dec;
    vNumCum *= dec;
    vDenCum *= dec;
  }

  void
  update(void)
  {
    valCum += val;
    vNumCum += vNum;
    vDenCum += vDen;
  }
};
}  // namespace incremental
}  // namespace EventEMin

#endif  // EVENT_EMIN_INCREMENTAL_DISPERSION_H
