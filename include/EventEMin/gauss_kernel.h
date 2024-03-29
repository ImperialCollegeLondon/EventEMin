#ifndef EVENT_EMIN_GAUSS_KERNEL_H
#define EVENT_EMIN_GAUSS_KERNEL_H

#include <cassert>
#include <cmath>

#include "EventEMin/types_def.h"

namespace EventEMin
{
namespace kernel
{
template <typename T>
T
gauss(const T& x)
{
  const T val = x * x;
  return std::exp(-T(0.5) * val);
}
template <typename T>
T
gauss(const T& x, const T& sigmaInv)
{
  const T val = x * sigmaInv * x;
  return std::exp(-T(0.5) * val);
}
template <typename T, int N>
T
gauss(const Ref<const Vector<T, N> >& x)
{
  const T val = x.transpose() * x;
  return std::exp(-T(0.5) * val);
}
template <typename T, int N>
T
gauss(const Ref<const Vector<T, N> >& x,
      const Ref<const Vector<T, N> >& sigmaInv)
{
  const T val = x.transpose() * (sigmaInv.asDiagonal() * x);
  return std::exp(-T(0.5) * val);
}
template <typename T>
T
gauss(const Ref<const Vector<T> >& x)
{
  const T val = x.transpose() * x;
  return std::exp(-T(0.5) * val);
}
template <typename T>
T
gauss(const Ref<const Vector<T> >& x, const Ref<const Vector<T> >& sigmaInv)
{
  assert(x.size() == sigmaInv.size());
  const T val = x.transpose() * (sigmaInv.asDiagonal() * x);
  return std::exp(-T(0.5) * val);
}

template <typename T, int N>
void
gaussKernelIterate(const int ksize, const int hksize, const int d, Vector<T>& p,
                   Array<Index, N>& ind, Tensor<T, N>& kernel)
{
  if (d >= N)
  {
    kernel(ind) = gauss<T>(p);
    return;
  }

  for (ind[d] = 0; ind[d] < ksize; ++ind[d])
  {
    p(d) = ind[d] - hksize;
    gaussKernelIterate<T, N>(ksize, hksize, d + 1, p, ind, kernel);
  }
}
template <typename T, int N>
void
gaussKernel(const int ksize, Tensor<T, N>& kernel)
{
  Array<Index, N> sizeArray;
  sizeArray.fill(ksize);
  kernel.resize(sizeArray);
  const int hksize = ksize >> 1;
  Vector<T> p(N);
  Array<Index, N> ind;
  gaussKernelIterate<T, N>(ksize, hksize, 0, p, ind, kernel);
  sizeArray.fill(hksize);
  kernel = kernel / (kernel(sizeArray) * std::pow(static_cast<T>(2.0 * M_PI),
                                                  static_cast<T>(N * 0.5)));
}
template <typename T>
void
gaussKernel(const int ksize, Vector<T>& kernel)
{
  kernel.resize(ksize);
  const int hksize = ksize >> 1;
  T p;
  for (int i = 0; i < ksize; ++i)
  {
    p = i - hksize;
    kernel(i) = gauss<T>(p);
  }
  kernel /= kernel(hksize, hksize) * std::sqrt(static_cast<T>(2.0 * M_PI));
}
template <typename T>
void
gaussKernel(const int ksize, Matrix<T>& kernel)
{
  kernel.resize(ksize, ksize);
  const int hksize = ksize >> 1;
  Vector<T, 2> p;
  for (int i = 0; i < ksize; ++i)
  {
    p(0) = i - hksize;
    for (int j = 0; j < ksize; ++j)
    {
      p(1) = j - hksize;
      kernel(i, j) = gauss<T, 2>(p);
    }
  }
  kernel /= kernel(hksize, hksize) * static_cast<T>(2.0 * M_PI);
}

template <typename T, int N>
void
gaussKernelIterate(const int ksize, const int hksize,
                   const Ref<const Vector<T> >& sigmaInv, const int d,
                   Vector<T>& p, Array<Index, N>& ind, Tensor<T, N>& kernel)
{
  if (d >= N)
  {
    kernel(ind) = gauss<T>(p, sigmaInv);
    return;
  }

  for (ind[d] = 0; ind[d] < ksize; ++ind[d])
  {
    p(d) = ind[d] - hksize;
    gaussKernelIterate<T, N>(ksize, hksize, sigmaInv, d + 1, p, ind, kernel);
  }
}
template <typename T, int N>
void
gaussKernel(const int ksize, const Ref<const Vector<T> >& sigmaInv,
            const T& sigmaDet, Tensor<T, N>& kernel)
{
  assert(sigmaInv.size() == N);

  Array<Index, N> sizeArray;
  sizeArray.fill(ksize);
  kernel.resize(sizeArray);
  const int hksize = ksize >> 1;
  Vector<T> p(N);
  Array<Index, N> ind;
  gaussKernelIterate<T, N>(ksize, hksize, sigmaInv, 0, p, ind, kernel);
  sizeArray.fill(hksize);
  kernel =
      kernel / (kernel(sizeArray) *
                std::pow(static_cast<T>(2.0 * M_PI), static_cast<T>(N * 0.5)) *
                std::sqrt(sigmaDet));
}
template <typename T>
void
gaussKernel(const int ksize, const T& sigmaInv, const T& sigmaDet,
            Vector<T>& kernel)
{
  kernel.resize(ksize);
  const int hksize = ksize >> 1;
  T p;
  for (int i = 0; i < ksize; ++i)
  {
    p = i - hksize;
    kernel(i) = gauss<T>(p, sigmaInv);
  }
  kernel /=
      kernel(hksize, hksize) * static_cast<T>(2.0 * M_PI) * std::sqrt(sigmaDet);
}
template <typename T>
void
gaussKernel(const int ksize, const Ref<const Vector<T> >& sigmaInv,
            const T& sigmaDet, Matrix<T>& kernel)
{
  assert(sigmaInv.size() == 2);

  kernel.resize(ksize, ksize);
  const int hksize = ksize >> 1;
  Vector<T, 2> p;
  for (int i = 0; i < ksize; ++i)
  {
    p(0) = i - hksize;
    for (int j = 0; j < ksize; ++j)
    {
      p(1) = j - hksize;
      kernel(i, j) = gauss<T, 2>(p, sigmaInv);
    }
  }
  kernel /=
      kernel(hksize, hksize) * static_cast<T>(2.0 * M_PI) * std::sqrt(sigmaDet);
}
}  // namespace kernel
}  // namespace EventEMin

#endif  // EVENT_EMIN_GAUSS_KERNEL_H
