#ifndef EVENT_EMIN_INCREMENTAL_POTENTIAL_H
#define EVENT_EMIN_INCREMENTAL_POTENTIAL_H

#include "EventEMin/dispersion/incremental_dispersion/dispersion.h"

namespace EventEMin
{
namespace incremental
{
template <typename M>
class Potential
{
 public:
  typedef M Model;
  typedef typename Model::T T;

  enum
  {
    NVars = Model::NVars,
    NDims = Model::NDims
  };

  typedef typename Model::Point Point;
  typedef typename Model::Vars Vars;
  typedef typename Model::CMatrix CMatrix;
  typedef typename Model::DMatrix DMatrix;
  typedef typename Model::PMatrix PMatrix;
  typedef Matrix<T, NVars, NDims> WMatrix;

 public:
  Potential(void) = default;

  void
  operator()(const Ref<const Point>& cmDiff, const Ref<const DMatrix>& dcmDiff,
             const Ref<const Point>& cmgDiff, const Ref<const PMatrix>& perDiff,
             VarsEstimate<Model>& varse) const
  {
    const T cmDiffExp = computeExp(T(-0.5) * (cmDiff.transpose() * cmDiff)(0));
    varse.val += cmDiffExp;
    const WMatrix cWeight(cmDiffExp * dcmDiff.transpose());
    varse.vNum.noalias() -= cWeight * cmgDiff;
    varse.vDen.noalias() += cWeight * perDiff;
  }
};
}  // namespace incremental

template <typename M>
using IncrementalPotential =
    incremental::Dispersion<incremental::Potential<M>,
                            incremental::Whitening<M, false> >;
template <typename M>
using IncrementalPotentialWhiten =
    incremental::Dispersion<incremental::Potential<M>,
                            incremental::Whitening<M, true> >;
}  // namespace EventEMin

#endif  // EVENT_EMIN_INCREMENTAL_POTENTIAL_H
