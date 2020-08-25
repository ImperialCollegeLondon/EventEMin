#ifndef ENTROPY_H
#define ENTROPY_H

#ifdef _OPENMP
  #include <csignal>
  #include <omp.h>
  #include <thread>
#endif

#include "dispersion/dispersion.h"

namespace event_model
{

namespace dispersion
{

template<typename Derived>
class Entropy : public Dispersion<Entropy<Derived> >
{

public:

  typedef typename DispersionTraits<Entropy<Derived> >::Model Model;
  typedef typename DispersionTraits<Entropy<Derived> >::T T;
  
private:

protected:

  const T cScale_;

public:

  Entropy(
    const T& cScale=T(100.0)):
    Dispersion<Entropy<Derived> >(),
    cScale_(cScale)
  {}

  T
  dimScale(
    const int d) const
  {
    return cScale_;
  }
  T
  halfOffset(void) const
  {
    return T(0.0);
  }
  T
  offset(void) const
  {
    return T(0.0);
  }

  template<typename U>
  U
  compute(
    const mtx<U>& c) const
  {
    U f=U(0.0);
#ifdef _OPENMP
    #pragma omp parallel shared(c, f)
#endif
    {
      const int npoints=this->npoints();
      mtx<U> cDiff;

#ifdef _OPENMP
    #pragma omp declare reduction(sum : U : omp_out+=omp_in) \
      initializer (omp_priv=U(0.0))
    #pragma omp for reduction(sum : f)
#endif
      for(int k=1; k<npoints; ++k)
      {
        pointsDiff(k, c, cDiff);
        f+=this->underlying().fCumulative(cDiff);
      }
    }

    return this->underlying().fFinal(f);
  }

protected:

  template<typename U>
  void
  pointsDiff(
    const int k,
    const mtx<U>& p,
    mtx<U>& pDiff) const
  {
    pDiff=p.leftCols(k).colwise()-p.col(k);
  }

private:

  Derived&
  underlying(void)
  {
    return static_cast<Derived&>(*this);
  }
  const Derived&
  underlying(void) const
  {
    return static_cast<const Derived&>(*this);
  }

};

}

}

#endif // ENTROPY_H
