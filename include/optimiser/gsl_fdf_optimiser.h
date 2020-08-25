#ifndef GSL_FDF_OPTIMIZER_H
#define GSL_FDF_OPTIMIZER_H

#include <cassert>
#include <Eigen/Core>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <unsupported/Eigen/AutoDiff>

#include "types_def.h"

namespace event_model
{

namespace optimiser
{

template<typename Func>
double
gslfdffFunc(
  const gsl_vector* gvars, void* params);

template<typename Func>
void
gslfdfdFunc(
  const gsl_vector* gvars, void* params,
  gsl_vector* gdf);

template<typename Func>
void
gslfdffdFunc(
  const gsl_vector* gvars, void* params,
  double* gf, gsl_vector* gdf);

struct GSLfdfOptimiserParams
{

  const gsl_multimin_fdfminimizer_type* type;
  const double iniStep, tol;
  const int maxIter, maxrIter, verbose;

  GSLfdfOptimiserParams(
    const gsl_multimin_fdfminimizer_type* type=gsl_multimin_fdfminimizer_conjugate_fr,
    const double iniStep=1.0e-1, const double tol=1.0e-6,
    const int maxIter=100, const int maxrIter=1, const int verbose=0):
    type(type),
    iniStep(iniStep), tol(tol),
    maxIter(maxIter), maxrIter(maxrIter), verbose(verbose)
  {}

};

template<typename F>
class GSLfdfOptimiser
{

public:

  typedef F Func;
  typedef typename Func::T T;
  typedef GSLfdfOptimiserParams OptimiserParams;

private:

  const OptimiserParams params_;

protected:

  const Eigen::AutoDiffJacobian<Func> func_;

  gsl_vector* gvars_;
  gsl_multimin_fdfminimizer* s_;
  gsl_multimin_function_fdf gslFunc_;

  Eigen::Map<vec<double> > gvarsMap_;
  vec<T, Func::nvars()> vars_;

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  GSLfdfOptimiser(
    const Func& func,
    const OptimiserParams& params=OptimiserParams()):
    params_(params),
    func_(func),
    gvars_(gsl_vector_alloc(nvars())),
    s_(gsl_multimin_fdfminimizer_alloc(params_.type, nvars())),
    gvarsMap_(gvars_->data, nvars())
  {
    gslFunc_.n=nvars();
    gslFunc_.f=&gslfdffFunc<Func>;
    gslFunc_.df=&gslfdfdFunc<Func>;
    gslFunc_.fdf=&gslfdffdFunc<Func>;
    gslFunc_.params=(void*) this;
  }
  GSLfdfOptimiser(
    const GSLfdfOptimiser& optimiser):
    params_(optimiser.params_),
    func_(optimiser.func_),
    gvars_(gsl_vector_alloc(nvars())),
    s_(gsl_multimin_fdfminimizer_alloc(params_.type, nvars())),
    gslFunc_(optimiser.gslFunc_),
    gvarsMap_(gvars_->data, nvars())
  {
    gslFunc_.params=(void*) this;
  }
  ~GSLfdfOptimiser(void)
  {
    if(gvars_!=nullptr)
    {
      gsl_vector_free(gvars_);
    }
    if(s_!=nullptr)
    {
      gsl_multimin_fdfminimizer_free(s_);
    }
  }

  static constexpr int
  nvars(void)
  {
    return Func::nvars();
  }
  const vec<T, nvars()>&
  vars(void) const
  {
    return vars_;
  }
  double
  iniStep(void) const
  {
    return params_.iniStep;
  }
  double
  tol(void) const
  {
    return params_.tol;
  }
  int
  maxIter(void) const
  {
    return params_.maxIter;
  }
  int
  maxrIter(void) const
  {
    return params_.maxrIter;
  }
  int
  verbose(void) const
  {
    return params_.verbose;
  }

  T
  run(
    const Eigen::Ref<const vec<T> >& iniVars)
  {
    assert(nvars()==iniVars.size());
    gvarsMap_=iniVars.template cast<double>();

    int iter=0, riter=0;
    int status;

    gsl_multimin_fdfminimizer_set(s_, &gslFunc_, gvars_, iniStep(), 1.0e-1);
    const Eigen::Map<const vec<double> > gvarsMap(s_->x->data, nvars());
    const Eigen::Map<const vec<double> > gradientMap(s_->gradient->data, nvars());

    do
    {
      ++iter;
      status=gsl_multimin_fdfminimizer_iterate(s_);
      if(status)
      {
        gsl_multimin_fdfminimizer_set(s_, &gslFunc_, s_->x, iniStep(), 1.0e-2);
        ++riter;
      }
      else
      {
        riter=0;
      }
      if(verbose())
      {
        std::cout<<iter<<", "<<riter<<": "<<s_->f<<' '<<gradientMap.norm()<<'\n';
        std::cout<<gradientMap.transpose()<<'\n';
      }
      status=gsl_multimin_test_gradient(s_->gradient, tol());
      if(verbose())
      {
        if(status==GSL_SUCCESS)
        {
          std::cout<<"minimum found at:\n";
        }
        std::cout<<gvarsMap.transpose()<<'\n';
      }
    } while(status==GSL_CONTINUE && iter<maxIter() && riter<maxrIter());

    vars_=gvarsMap.cast<T>();
    return static_cast<T>(s_->f);
  }

  double
  operator()(
    const gsl_vector* gvars) const
  {
    const Eigen::Map<const vec<double> > varsMap(gvars->data, nvars());
    const vec<T, nvars()> vars(varsMap.cast<T>());
    vec<T, 1> f;
    func_(vars, &f, nullptr);
    return static_cast<double>(f(0));
  }

  void
  operator()(
    const gsl_vector* gvars,
    gsl_vector* gdf) const
  {
    const Eigen::Map<const vec<double> > varsMap(gvars->data, nvars());
    Eigen::Map<vec<double> > dfMap(gdf->data, nvars());
    const vec<T, nvars()> vars(varsMap.cast<T>());
    vec<T, 1> f;
    mtx<T, 1, nvars()> df;
    func_(vars, &f, &df);
    dfMap=df.template cast<double>();
  }

  void
  operator()(
    const gsl_vector* gvars,
    double* gf, gsl_vector* gdf) const
  {
    const Eigen::Map<const vec<double> > varsMap(gvars->data, nvars());
    Eigen::Map<vec<double> > dfMap(gdf->data, nvars());
    const vec<T, nvars()> vars(varsMap.cast<T>());
    vec<T, 1> f;
    mtx<T, 1, nvars()> df;
    func_(vars, &f, &df);
    *gf=static_cast<double>(f(0));
    dfMap=df.template cast<double>();
  }

protected:

private:
  
};

template<typename Func>
double
gslfdffFunc(
  const gsl_vector* gvars, void* params)
{
  return (*static_cast<GSLfdfOptimiser<Func>*>(params))(gvars);
}

template<typename Func>
void
gslfdfdFunc(
  const gsl_vector* gvars, void* params,
  gsl_vector* gdf)
{
  (*static_cast<GSLfdfOptimiser<Func>*>(params))(gvars, gdf);
}

template<typename Func>
void
gslfdffdFunc(
  const gsl_vector* gvars, void* params,
  double* gf, gsl_vector* gdf)
{
  (*static_cast<GSLfdfOptimiser<Func>*>(params))(gvars, gf, gdf);
}

}

}

#endif // GSL_OPTIMIZER_H
