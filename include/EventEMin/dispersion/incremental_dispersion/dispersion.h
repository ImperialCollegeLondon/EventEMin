#ifndef EVENT_EMIN_INCREMENTAL_DISPERSION_IMPL_H
#define EVENT_EMIN_INCREMENTAL_DISPERSION_IMPL_H

#include "EventEMin/dispersion/incremental_dispersion.h"

namespace EventEMin
{
namespace incremental
{
template <typename DispersionImpl, typename Whitening>
class Dispersion
{
 public:
  typedef typename DispersionImpl::Model Model;
  typedef typename DispersionImpl::T T;

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

  typedef DispersionParams<T> Params;

  typedef StdVector<Cell> Grid;

 private:
  Matrix<T, 3, 3> camParams_;
  projectEvent<T, NDims> projectEvent_;

  Point scale_;

  Array<Index, 2> cMin_, cMax_;
  Array<Index, 2> dim_;

  Params params_;

  int nPointsCur_, nPointsMax_, nBuffer_;
  int curInd_, prevInd_, refInd_, nPoints_;

  StdVector<int> inds_;
  StdVector<Point> c_;
  StdVector<T> ts_;
  StdVector<T> tsDiffRef_;
  StdVector<Vector<int, 2> > ij_;

  T nDecay_;

  int iter_;

  DispersionImpl dispersionImpl_;
  Whitening whitening_;

 protected:
  Grid grid_;

  Model model_;
  VarsEstimate<Model> varse_;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Dispersion(void) = delete;
  Dispersion(const Matrix<T, 3, 3>& camParams, const Point& scale,
             const Params& params, const int nPoints,
             const Array<Index, 2>& dim)
      : camParams_(camParams),
        scale_(scale),
        dim_(dim),
        params_(params),
        nPointsCur_(nPoints),
        nPointsMax_(nPoints),
        nBuffer_(nextPower2(nPointsMax_)),
        curInd_(0),
        prevInd_(0),
        refInd_(0),
        nPoints_(0),
        c_(nBuffer_),
        ts_(nBuffer_),
        ij_(nBuffer_),
        nDecay_(T(0.0)),
        iter_(0)
  {
    --nBuffer_;
    int gridSize = 1;
    for (int d = 0; d < 2; ++d)
    {
      cMin_[d] = 0;
      cMax_[d] = dim_[d] - 1;
      gridSize *= dim_[d];
    }
    grid_.resize(gridSize, Cell());
  }

  T
  minStep(void) const
  {
    return params_.minStep;
  }
  int
  maxIter(void) const
  {
    return params_.maxIter;
  }
  int
  wSize(void) const
  {
    return params_.wSize;
  }
  int
  iter(void) const
  {
    return iter_;
  }
  const Vars&
  vars(void) const
  {
    return varse_.vars;
  }

  void
  decay(const T& ts)
  {
    const T dec = computeExp(-nDecay_ * (ts - ts_[prevInd_]));
    varse_.decay(dec);
    nDecay_ *= dec;
  }

  void
  run(const Ref<const Point>& c, const T& ts)
  {
    Point cu;
    projectEvent_(camParams_, c, cu);
    const Vector<int, 2> ij(
        cu.template head<2>().array().round().template cast<int>());
    const int ind = ij2ind(ij);
    Vector<int, 2> ijMin, ijMax;
    ijBoundary(ij, ijMin, ijMax);
    ij2ind(ijMin, ijMax, grid_, inds_);

    // circular shift
    T nProp;
    if (nPoints_ < nPointsCur_)
    {
      ++nPoints_;
      nProp = T(nPoints_) / nPointsCur_;
    }
    else
    {
      nProp = T(1.0);
      refInd_ = fastIncrementCyclic(refInd_, nBuffer_);
      // remove old event
      grid_[ij2ind(ij_[refInd_])].remove();
    }
    // compute difference in times
    tsDiffRef_.clear();
    for (int i = 0; i < static_cast<int>(inds_.size()); ++i)
    {
      tsDiffRef_.push_back(ts_[inds_[i]] - ts_[refInd_]);
    }
    const T tsDiffRef = ts - ts_[refInd_];

    // decay factor
    decay(ts);

    // compute the estimate with the new event
    T d = minStep();
    for (iter_ = 0; iter_ < maxIter() && d >= minStep(); ++iter_)
    {
      // std::cout << "iter: " << iter_ << '\n';
      d = iterate(inds_, nProp, c, tsDiffRef, varse_);
    }
    // update variables
    varse_.update();
    ++nDecay_;

    // add event to grid
    grid_[ind].add(curInd_);
    c_[curInd_] = c;
    ts_[curInd_] = ts;
    ij_[curInd_] = ij;

    // circular shift
    prevInd_ = curInd_;
    curInd_ = fastIncrementCyclic(curInd_, nBuffer_);

    // update stats incrementally
    whitening_.updateStats(varse_.vars, c, tsDiffRef, nPoints_);
  }

  void
  setInc(const T& inc)
  {
    const int nPointsCur = T(nPointsMax_) / inc;
    if (nPoints_ >= nPointsCur)
    {
      nPoints_ = nPointsCur;
      for (int indDiff = nPointsCur_ - nPointsCur; indDiff > 0; --indDiff)
      {
        refInd_ = fastIncrementCyclic(refInd_, nBuffer_);
        // remove old event
        grid_[ij2ind(ij_[refInd_])].remove();
      }
    }
    nPointsCur_ = nPointsCur;
  }

 protected:
  T
  iterate(const StdVector<int>& inds, const T& nProp, const Ref<const Point>& c,
          const T& tsDiffRef, VarsEstimate<Model>& varse) const
  {
    Point cm, cml, cmDiff, cg, cgl, cms, cmsl;
    DMatrix dcm, dcml;
    PMatrix per, perl;
    model_(varse.vars, c, tsDiffRef, cm, dcm, cg, per);
    whitening_.processEvent(scale_, cm, cms);

    varse.val = T(0.0);
    varse.vNum.setZero();
    varse.vDen.setZero();

    for (int i = 0; i < static_cast<int>(inds.size()); ++i)
    {
      model_(varse.vars, c_[inds[i]], tsDiffRef_[i], cml, dcml, cgl, perl);
      whitening_.processEvent(scale_, cml, cmsl);
      cmDiff = cms - cmsl;

      dispersionImpl_(cmDiff, dcm - dcml, cm - cml - cg + cgl, per - perl,
                      varse);
    }

    if (T(0.0) < varse.val)
    {
      varse.vNum *= nProp;  // for stability with few points
      // least-squares solution
      Vars v((varse.vDenCum + varse.vDen)
                 .template selfadjointView<Eigen::Upper>()
                 .ldlt()  // pos/neg semidefinite
                 //.partialPivLu()  // invertible
                 //.fullPivLu()
                 .solve(varse.vNumCum + varse.vNum));
      v -= varse.vars;
      varse.vars += v;
      const T d = v.norm();

      /*std::cout << "val: " << varse.val << ' ' << nProp << '\n';
      std::cout << "vNum: " << varse.vNum.transpose() << '\n';
      std::cout << "vNumCum: " << (varse.vNumCum + varse.vNum).transpose()
                << '\n';
      std::cout << "vDen: " << varse.vDen << '\n';
      std::cout << "vDenCum: " << varse.vDenCum + varse.vDen << '\n';
      std::cout << "vars: " << varse.vars.transpose() << '\n';
      std::cout << "v: " << v.transpose() << ", d: " << d << '\n';*/

      return d;
    }
    return T(0.0);
  }

 private:
  int
  ij2ind(const Ref<const Vector<int, 2> >& ij) const
  {
    return ij(0) * dim_[1] + ij(1);
  }
  void
  ij2ind(const Ref<const Vector<int, 2> >& ijMin,
         const Ref<const Vector<int, 2> >& ijMax, const Grid& grid,
         StdVector<int>& inds) const
  {
    inds.clear();
    Vector<int, 2> ij;
    for (ij(0) = ijMin(0); ij(0) <= ijMax(0); ++ij(0))
    {
      for (ij(1) = ijMin(1); ij(1) <= ijMax(1); ++ij(1))
      {
        const int ind = ij2ind(ij);
        for (int i = 0; i < grid[ind].n(); ++i)
        {
          inds.push_back(grid[ind].ind(i));
        }
      }
    }
  }

  void
  ijBoundary(const Ref<const Vector<int, 2> >& ij, Vector<int, 2>& ijMin,
             Vector<int, 2>& ijMax) const
  {
    for (int d = 0; d < 2; ++d)
    {
      ijMin(d) = std::max(ij(d) - wSize(), static_cast<int>(cMin_[d]));
      ijMax(d) = std::min(ij(d) + wSize(), static_cast<int>(cMax_[d]));
    }
  }
};
}  // namespace incremental
}  // namespace EventEMin

#endif  // EVENT_EMIN_INCREMENTAL_DISPERSION_IMPL_H
