#ifndef EVENT_TRANSFORM_H
#define EVENT_TRANSFORM_H

#include <cassert>
#include <cmath>
#include <Eigen/Core>

#include "types_def.h"

namespace event
{

template<typename T>
void
computeBoundaries(
  const Eigen::Ref<const vec<T> >& c,
  const Eigen::Ref<const vec<int> >& cMin, const Eigen::Ref<const vec<int> >& cMax,
  vec<int>& lMin, vec<int>& lMax)
{
  assert(c.size()==cMin.size() && c.size()==cMax.size());
  lMin.resize(c.size()), lMax.resize(c.size());
  for(int d=0; d<c.size(); ++d)
  {
    lMin(d)=std::max(static_cast<int>(std::floor(c(d))), cMin(d));
    lMax(d)=std::min(static_cast<int>(std::ceil(c(d))), cMax(d));
  }
}

template<typename T>
void
projectEvent(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  const Eigen::Ref<const vec<T> >& c3,
  Eigen::Ref<vec<T> > c2)
{
  assert(c3.size()==3 && c2.size()==2);
  c2(0)<<camParams(0, 2)+camParams(0, 0)*c3(0)/c3(2),
    camParams(1, 2)-camParams(1, 1)*c3(1)/c3(2);
}
template<typename T>
void
projectEvents(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  const Eigen::Ref<const mtx<T> >& c3,
  mtx<T>& c2)
{
  assert(c3.rows()==3 && c2.rows()==2 && c3.cols()==c2.cols());
  c2.row(0)=camParams(0, 2)+camParams(0, 0)*c3.row(0).array()/c3.row(2).array();
  c2.row(1)=camParams(1, 2)-camParams(1, 1)*c3.row(1).array()/c3.row(2).array();
}

template<typename T>
void
unprojectEvent(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  const Eigen::Ref<const vec<T> >& c2,
  const T& depthInv,
  Eigen::Ref<vec<T> > c3)
{
  assert(c2.size()==2 && c3.size()==3);
  c3<<(c2(0)-camParams(0, 2))/(camParams(0, 0)*depthInv),
    (camParams(1, 2)-c2(1))/(camParams(1, 1)*depthInv),
    T(1.0)/depthInv;
}
template<typename T>
void
unprojectEvents(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  const Eigen::Ref<const mtx<T> >& c2,
  const Eigen::Ref<const vec<T> >& depthInv,
  mtx<T>& c3)
{
  assert(c2.rows()==2 && c3.rows()==3);
  assert(depthInv.size()==c2.cols() && c3.cols()==c2.cols());
  c3.row(0)=(c2.row(0).array()-camParams(0, 2))/(camParams(0, 0)*depthInv.transpose().array());
  c3.row(1)=(camParams(1, 2)-c2.row(1).array())/(camParams(1, 1)*depthInv.transpose().array());
  c3.row(2)=depthInv.transpose().array().inverse();
}

template<typename T>
void
unprojectEvent(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  Eigen::Ref<vec<T> > c)
{
  assert(c.size()==3);
  c<<c(2)*(c(0)-camParams(0, 2))/camParams(0, 0),
    c(2)*(camParams(1, 2)-c(1))/camParams(1, 1);
}
template<typename T>
void
unprojectEvents(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  mtx<T>& c)
{
  assert(c.rows()==3);
  c.row(0)=c.row(2).array()*(c.row(0).array()-camParams(0, 2))/camParams(0, 0);
  c.row(1)=c.row(2).array()*(camParams(1, 2)-c.row(1).array())/camParams(1, 1);
}

template<typename T>
void
transformEvent(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  Eigen::Ref<vec<T> > c)
{
  assert(c.size()==2);
  c<<(c(0)-camParams(0, 2))/camParams(0, 0),
    (camParams(1, 2)-c(1))/camParams(1, 1);
}
template<typename T>
void
transformEvents(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  mtx<T>& c)
{
  assert(c.rows()==2);
  c.row(0)=(c.row(0).array()-camParams(0, 2))/camParams(0, 0);
  c.row(1)=(camParams(1, 2)-c.row(1).array())/camParams(1, 1);
}

template<typename T>
void
transformEvent(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  const Eigen::Ref<const vec<T> >& c,
  Eigen::Ref<vec<T> > ct)
{
  assert(c.size()==2 && ct.size()==2);
  ct<<(c(0)-camParams(0, 2))/camParams(0, 0),
    (camParams(1, 2)-c(1))/camParams(1, 1);
}
template<typename T>
void
transformEvents(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  const Eigen::Ref<const mtx<T> >& c,
  mtx<T>& ct)
{
  assert(c.rows()==2 && ct.rows()==2 && c.cols()==ct.cols());
  ct.row(0)=(c.row(0).array()-camParams(0, 2))/camParams(0, 0);
  ct.row(1)=(camParams(1, 2)-c.row(1).array())/camParams(1, 1);
}

template<typename T>
void
untransformEvent(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  Eigen::Ref<vec<T> > ct)
{
  assert(ct.size()==2);
  ct<<camParams(0, 2)+camParams(0, 0)*ct(0),
    camParams(1, 2)-camParams(1, 1)*ct(1);
}
template<typename T>
void
untransformEvents(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  mtx<T>& ct)
{
  assert(ct.rows()==2);
  ct.row(0)=camParams(0, 2)+camParams(0, 0)*ct.row(0).array();
  ct.row(1)=camParams(1, 2)-camParams(1, 1)*ct.row(1).array();
}

template<typename T>
void
untransformEvent(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  const Eigen::Ref<const vec<T> >& ct,
  Eigen::Ref<vec<T> > c)
{
  assert(c.size()==2 && ct.size()==2);
  c<<camParams(0, 2)+camParams(0, 0)*ct(0),
    camParams(1, 2)-camParams(1, 1)*ct(1);
}
template<typename T>
void
untransformEvents(
  const Eigen::Ref<const mtx<T, 3, 3> >& camParams,
  const Eigen::Ref<const mtx<T> >& ct,
  mtx<T>& c)
{
  assert(c.rows()==2 && ct.rows()==2 && c.cols()==ct.cols());
  c.row(0)=camParams(0, 2)+camParams(0, 0)*ct.row(0).array();
  c.row(1)=camParams(1, 2)-camParams(1, 1)*ct.row(1).array();
}

}

#endif // EVENT_TRANSFORM_H
