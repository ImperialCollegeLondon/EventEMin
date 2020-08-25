#ifndef GAUSS_KERNEL_H
#define GAUSS_KERNEL_H

#include <cassert>
#include <cmath>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

#include "types_def.h"

namespace kernel
{

template<typename T>
T
gauss(
  const T& x)
{
  const T val=x*x;
  return std::exp(-T(0.5)*val);
}
template<typename T>
T
gauss(
  const T& x,
  const T& sigmaInv)
{
  const T val=x*sigmaInv*x;
  return std::exp(-T(0.5)*val);
}
template<typename T, int N>
T
gauss(
  const Eigen::Ref<const vec<T, N> >& x)
{
  const T val=x.transpose()*x;
  return std::exp(-T(0.5)*val);
}
template<typename T, int N>
T
gauss(
  const Eigen::Ref<const vec<T, N> >& x,
  const Eigen::Ref<const vec<T, N> >& sigmaInv)
{
  const T val=x.transpose()*(sigmaInv.asDiagonal()*x);
  return std::exp(-T(0.5)*val);
}
template<typename T>
T
gauss(
  const Eigen::Ref<const vec<T> >& x)
{
  const T val=x.transpose()*x;
  return std::exp(-T(0.5)*val);
}
template<typename T>
T
gauss(
  const Eigen::Ref<const vec<T> >& x,
  const Eigen::Ref<const vec<T> >& sigmaInv)
{
  assert(x.size()==sigmaInv.size());
  const T val=x.transpose()*(sigmaInv.asDiagonal()*x);
  return std::exp(-T(0.5)*val);
}

template<typename T>
T
dgauss(
  const T& g,
  const T& x)
{
  return -g*x;
}
template<typename T>
T
dgauss(
  const T& g,
  const T& x,
  const T& sigmaInv)
{
  return -g*x*sigmaInv;
}
template<typename T, int N>
void
dgauss(
  const T& g,
  const Eigen::Ref<const vec<T, N> >& x,
  Eigen::Ref<vec<T, N> > dg)
{
  dg=-g*x;
}
template<typename T, int N>
void
dgauss(
  const T& g,
  const Eigen::Ref<const vec<T, N> >& x,
  const Eigen::Ref<const vec<T, N> >& sigmaInv,
  Eigen::Ref<vec<T, N> > dg)
{
  dg=-g*(x.transpose()*sigmaInv.asDiagonal()).transpose();
}
template<typename T>
void
dgauss(
  const T& g,
  const Eigen::Ref<const vec<T> >& x,
  Eigen::Ref<vec<T> > dg)
{
  dg=-g*x;
}
template<typename T>
void
dgauss(
  const T& g,
  const Eigen::Ref<const vec<T> >& x,
  const Eigen::Ref<const vec<T> >& sigmaInv,
  Eigen::Ref<vec<T> > dg)
{
  assert(x.size()==sigmaInv.size());
  dg=-g*(x.transpose()*sigmaInv.asDiagonal()).transpose();
}

template<typename T, int NDims>
void
gaussKernelIterate(
  const int ksize, const int hksize,
  const int d,
  vec<T>& p,
  array<int, NDims>& ind,
  tensor<T, NDims>& kernel)
{
  if(d>=NDims)
  {
    kernel(ind)=gauss<T>(p);
    return;
  }

  for(ind.at(d)=0; ind.at(d)<ksize; ++ind.at(d))
  {
    p(d)=ind.at(d)-hksize;
    gaussKernelIterate<T, NDims>(ksize, hksize, d+1, p, ind, kernel);
  }
}
template<typename T, int NDims>
void
gaussKernel(
  const int ksize,
  tensor<T, NDims>& kernel)
{
  array<int, NDims> sizeArray;
  sizeArray.fill(ksize);
  kernel.resize(sizeArray);
  const int hksize=ksize>>1;
  vec<T> p(NDims);
  array<int, NDims> ind;
  gaussKernelIterate<T, NDims>(ksize, hksize, 0, p, ind, kernel);
  sizeArray.fill(hksize);
  kernel=kernel/(kernel(sizeArray)*static_cast<T>(std::pow(2.0*M_PI, static_cast<T>(NDims)*0.5)));
}
template<typename T>
void
gaussKernel(
  const int ksize,
  vec<T>& kernel)
{
  kernel.resize(ksize);
  const int hksize=ksize>>1;
  T p;
  for(int i=0; i<ksize; ++i)
  {
    p=i-hksize;
    kernel(i)=gauss<T>(p);
  }
  kernel/=kernel(hksize, hksize)*std::sqrt(2.0*M_PI);
}
template<typename T>
void
gaussKernel(
  const int ksize,
  mtx<T>& kernel)
{
  kernel.resize(ksize, ksize);
  const int hksize=ksize>>1;
  vec<T, 2> p;
  for(int i=0; i<ksize; ++i)
  {
    p(0)=i-hksize;
    for(int j=0; j<ksize; ++j)
    {
      p(1)=j-hksize;
      kernel(i, j)=gauss<T, 2>(p);
    }
  }
  kernel/=kernel(hksize, hksize)*2.0*M_PI;
}

template<typename T, int NDims>
void
gaussKernelIterate(
  const int ksize, const int hksize,
  const Eigen::Ref<const vec<T> >& sigmaInv,
  const int d,
  vec<T>& p,
  array<int, NDims>& ind,
  tensor<T, NDims>& kernel)
{
  if(d>=NDims)
  {
    kernel(ind)=gauss<T>(p, sigmaInv);
    return;
  }

  for(ind.at(d)=0; ind.at(d)<ksize; ++ind.at(d))
  {
    p(d)=ind.at(d)-hksize;
    gaussKernelIterate<T, NDims>(ksize, hksize, sigmaInv, d+1, p, ind, kernel);
  }
}
template<typename T, int NDims>
void
gaussKernel(
  const int ksize,
  const Eigen::Ref<const vec<T> >& sigmaInv,
  const T& sigmaDet,
  tensor<T, NDims>& kernel)
{
  assert(sigmaInv.size()==NDims);

  array<int, NDims> sizeArray;
  sizeArray.fill(ksize);
  kernel.resize(sizeArray);
  const int hksize=ksize>>1;
  vec<T> p(NDims);
  array<int, NDims> ind;
  gaussKernelIterate<T, NDims>(ksize, hksize, sigmaInv, 0, p, ind, kernel);
  sizeArray.fill(hksize);
  kernel=kernel/(kernel(sizeArray)*static_cast<T>(std::pow(2.0*M_PI, static_cast<T>(NDims)*0.5)*std::sqrt(sigmaDet)));
}
template<typename T>
void
gaussKernel(
  const int ksize,
  const T& sigmaInv,
  const T& sigmaDet,
  vec<T>& kernel)
{
  kernel.resize(ksize);
  const int hksize=ksize>>1;
  T p;
  for(int i=0; i<ksize; ++i)
  {
    p=i-hksize;
    kernel(i)=gauss<T>(p, sigmaInv);
  }
  kernel/=kernel(hksize, hksize)*2.0*M_PI*std::sqrt(sigmaDet);
}
template<typename T>
void
gaussKernel(
  const int ksize,
  const Eigen::Ref<const vec<T> >& sigmaInv,
  const T& sigmaDet,
  mtx<T>& kernel)
{
  assert(sigmaInv.size()==2);
  kernel.resize(ksize, ksize);
  const int hksize=ksize>>1;
  vec<T, 2> p;
  for(int i=0; i<ksize; ++i)
  {
    p(0)=i-hksize;
    for(int j=0; j<ksize; ++j)
    {
      p(1)=j-hksize;
      kernel(i, j)=gauss<T, 2>(p, sigmaInv);
    }
  }
  kernel/=kernel(hksize, hksize)*2.0*M_PI*std::sqrt(sigmaDet);
}

}

#endif // GAUSS_KERNEL_H
