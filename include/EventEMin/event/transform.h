#ifndef EVENT_EMIN_EVENT_TRANSFORM_H
#define EVENT_EMIN_EVENT_TRANSFORM_H

#include <cassert>
#include <cmath>

#include "EventEMin/types_def.h"

namespace EventEMin
{
// partial specialisation only works for structs not functions
template <typename I, int N, typename T>
struct computeBoundaries
{
  void
  operator()(const Ref<const Vector<T, N> >& c, const Array<I, N>& cMin,
             const Array<I, N>& cMax, Array<I, N>& lMin,
             Array<I, N>& lMax) const
  {
    for (I d = 0; d < N; ++d)
    {
      lMin[d] = std::max(static_cast<I>(std::floor(c(d).value())), cMin[d]);
      lMax[d] = std::min(static_cast<I>(std::ceil(c(d).value())), cMax[d]);
    }
  }
};
template <typename I, int N>
struct computeBoundaries<I, N, double>
{
  void
  operator()(const Ref<const Vector<double, N> >& c, const Array<I, N>& cMin,
             const Array<I, N>& cMax, Array<I, N>& lMin,
             Array<I, N>& lMax) const
  {
    for (I d = 0; d < N; ++d)
    {
      lMin[d] = std::max(static_cast<I>(std::floor(c(d))), cMin[d]);
      lMax[d] = std::min(static_cast<I>(std::ceil(c(d))), cMax[d]);
    }
  }
};
template <typename I, int N>
struct computeBoundaries<I, N, float>
{
  void
  operator()(const Ref<const Vector<float, N> >& c, const Array<I, N>& cMin,
             const Array<I, N>& cMax, Array<I, N>& lMin,
             Array<I, N>& lMax) const
  {
    for (I d = 0; d < N; ++d)
    {
      lMin[d] = std::max(static_cast<I>(std::floor(c(d))), cMin[d]);
      lMax[d] = std::min(static_cast<I>(std::ceil(c(d))), cMax[d]);
    }
  }
};

template <typename T, int N>
struct projectEvent;
template <typename T>
struct projectEvent<T, 2>
{
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             Ref<Vector<T> > c) const
  {
    assert(c.size() >= 2);

    c(0) = camParams(0, 0) * c(0) + camParams(0, 2);
    c(1) = camParams(1, 1) * c(1) + camParams(1, 2);
  }
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             const Ref<const Vector<T> >& cIn, Ref<Vector<T> > cOut) const
  {
    assert(cIn.size() >= 2);
    assert(cOut.size() >= 2);

    cOut(0) = camParams(0, 0) * cIn(0) + camParams(0, 2);
    cOut(1) = camParams(1, 1) * cIn(1) + camParams(1, 2);
  }
};
template <typename T>
struct projectEvent<T, 3>
{
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             Ref<Vector<T> > c) const
  {
    assert(c.size() >= 3);

    c.template head<3>() = camParams * c.template head<3>();
    c.template head<2>() /= c(2);
  }
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             const Ref<const Vector<T> >& cIn, Ref<Vector<T> > cOut) const
  {
    assert(cIn.size() >= 3);
    assert(cOut.size() >= 3);

    cOut.template head<3>().noalias() = camParams * cIn.template head<3>();
    cOut.template head<2>() /= cIn(2);
  }
};

template <typename T, int N>
struct projectEvents;
template <typename T>
struct projectEvents<T, 2>
{
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             Ref<Matrix<T> > c) const
  {
    assert(c.rows() >= 2);

    c.row(0) = camParams(0, 0) * c.row(0).array() + camParams(0, 2);
    c.row(1) = camParams(1, 1) * c.row(1).array() + camParams(1, 2);
  }
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             const Ref<const Matrix<T> >& cIn, Ref<Matrix<T> > cOut) const
  {
    assert(cIn.rows() >= 2);
    assert(cOut.rows() >= 2);
    assert(cIn.cols() == cOut.cols());

    cOut.row(0) = camParams(0, 0) * cIn.row(0).array() + camParams(0, 2);
    cOut.row(1) = camParams(1, 1) * cIn.row(1).array() + camParams(1, 2);
  }
};
template <typename T>
struct projectEvents<T, 3>
{
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             Ref<Matrix<T> > c) const
  {
    assert(c.rows() >= 3);

    c.template topRows<3>() = camParams * c.template topRows<3>();
    c.template topRows<2>().rowwise().array() /= c.row(2).array();
  }
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             const Ref<const Matrix<T> >& cIn, Ref<Matrix<T> > cOut) const
  {
    assert(cIn.rows() >= 3);
    assert(cOut.rows() >= 3);
    assert(cIn.cols() == cOut.cols());

    cOut.template topRows<3>().noalias() =
        camParams * cIn.template topRows<3>();
    cOut.template topRows<2>().array().rowwise() /= cIn.row(2).array();
  }
};

template <typename T, int N>
struct unprojectEvent;
template <typename T>
struct unprojectEvent<T, 2>
{
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             Ref<Vector<T> > c) const
  {
    assert(c.size() >= 2);

    c(0) = (c(0) - camParams(0, 2)) / camParams(0, 0);
    c(1) = (c(1) - camParams(1, 2)) / camParams(1, 1);
  }
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             const Ref<const Vector<T> >& cIn, Ref<Vector<T> > cOut) const
  {
    assert(cIn.size() >= 2);
    assert(cOut.size() >= 2);

    cOut(0) = (cIn(0) - camParams(0, 2)) / camParams(0, 0);
    cOut(1) = (cIn(1) - camParams(1, 2)) / camParams(1, 1);
  }
};
template <typename T>
struct unprojectEvent<T, 3>
{
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             Ref<Vector<T> > c) const
  {
    assert(c.size() >= 3);

    c.template head<2>() *= c(2);
    c.template head<3>() = camParams.inverse() * c.template head<3>();
  }
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             const Ref<const Vector<T> >& cIn, Ref<Vector<T> > cOut) const
  {
    assert(cIn.size() >= 3);
    assert(cOut.size() >= 3);

    cOut.template head<3>() = cIn.template head<3>();
    this->operator()(camParams, cOut);
  }
};

template <typename T, int N>
struct unprojectEvents;
template <typename T>
struct unprojectEvents<T, 2>
{
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             Ref<Matrix<T> > c) const
  {
    assert(c.rows() >= 2);

    c.row(0) = (c.row(0).array() - camParams(0, 2)) / camParams(0, 0);
    c.row(1) = (c.row(1).array() - camParams(1, 2)) / camParams(1, 1);
  }
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             const Ref<const Matrix<T> >& cIn, Ref<Matrix<T> > cOut) const
  {
    assert(cIn.rows() >= 2);
    assert(cOut.rows() >= 2);
    assert(cIn.cols() == cOut.cols());

    cOut.row(0) = (cIn.row(0).array() - camParams(0, 2)) / camParams(0, 0);
    cOut.row(1) = (cIn.row(1).array() - camParams(1, 2)) / camParams(1, 1);
  }
};
template <typename T>
struct unprojectEvents<T, 3>
{
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             Ref<Matrix<T> > c) const
  {
    assert(c.rows() >= 3);

    c.template topRows<2>().array().rowwise() *= c.row(2).array();
    c.template topRows<3>() = camParams.inverse() * c.template topRows<3>();
  }
  void
  operator()(const Ref<const Matrix<T, 3, 3> >& camParams,
             const Ref<const Matrix<T> >& cIn, Ref<Matrix<T> > cOut) const
  {
    assert(cIn.rows() >= 3);
    assert(cOut.rows() >= 3);
    assert(cIn.cols() == cOut.cols());

    cOut.template topRows<3>() = cIn.template topRows<3>();
    this->operator()(camParams, cOut);
  }
};
}  // namespace EventEMin

#endif  // EVENT_EMIN_EVENT_TRANSFORM_H
