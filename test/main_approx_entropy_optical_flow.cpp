#include <iostream>
#include <string>

#include "dispersion.h"
#include "event.h"
#include "image.h"
#include "model.h"
#include "optimiser.h"
#include "types_def.h"

using namespace event_model;

int
main(void)
{
  typedef float T;

  // model
  typedef model::OpticalFlow<T> Model;

  constexpr int ndims=Model::ndims(), nvars=Model::nvars();

  // you can modify the entropy used by uncommenting the corresponding line
  //typedef dispersion::ApproxPotential<Model> Dispersion;
  //typedef dispersion::ApproxRenyi<Model> Dispersion;
  //typedef dispersion::ApproxShannon<Model> Dispersion;
  //typedef dispersion::ApproxSharmaMittal<Model> Dispersion;
  typedef dispersion::ApproxTsallis<Model> Dispersion;

  // optimiser
  typedef optimiser::GSLfdfOptimiser<Dispersion> Optimiser;

  const std::string fdir("../dataset/poster_translation");

  // read events from file
  event::Events<T, ndims> evs;
  const std::string fevents(fdir+"/events.txt");
  if(!event::load<T, ndims>(fevents, evs))
  {
    std::cout<<"error reading events from file "<<fevents<<'\n';
    return -1;
  }
  const int nEvents=evs.size();
  std::cout<<"number of events: "<<nEvents<<'\n';

  int width, height;
  mtx<T, 3, 3> camParams;
  const std::string fcalib(fdir+"/calib.txt");
  if(!loadCamParams<T>(fcalib, width, height, camParams))
  {
    std::cout<<"error reading calibration from file "<<fcalib<<'\n';
    return -1;
  }

  mtx<T> c;
  vec<T> ts;
  vec<int> polarity;
  event::events2eigen<T>(evs, c, ts, polarity);

  // show original events
  showGray<T>(c, polarity, width, height, "Original Events");
  std::cout<<"press any key to continue...\n";
  cv::waitKey(0);

  mtx<T> ct(ndims, nEvents);
  event::transformEvents<T>(camParams, c, ct);

  Dispersion dispersion({width, height});
  dispersion.assignPoints(ct, ts, polarity);

  // initial parameters
  vec<T, nvars> vars;
  vars.setConstant(1.0e-6);

  // initial step size of the optimisation
  const double iniStep=1.0;
  // tolerance that indicates a minimum has been reached
  const double tol=1.0e-16;
  // maximum iterations
  const int maxIter=100, maxrIter=1;
  // optimiser status feedback: 0 - no feedback, 1 - feedback at each iteration
  const int verbosity=1;
  // optimise
  Optimiser optimiser(dispersion, Optimiser::OptimiserParams(
    gsl_multimin_fdfminimizer_conjugate_fr, iniStep, tol, maxIter, maxrIter, verbosity));
  std::cout<<"score: "<<optimiser.run(vars)<<", v: "<<optimiser.vars().transpose()<<'\n';

  // show modelled events according to the estimated parameters
  const Model model;
  mtx<T> cm(ndims, nEvents), ctm(ndims, nEvents);
  model(optimiser.vars(), ct, ts, ctm);
  event::untransformEvents<T>(camParams, ctm, cm);
  showGray<T>(cm, polarity, width, height, "Modelled Events");
  std::cout<<"press any key to exit...\n";
  cv::waitKey(0);

  return 0;
}
