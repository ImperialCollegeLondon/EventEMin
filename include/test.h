#include <iostream>
#include <string>

#include "event.h"
#include "image.h"
#include "optimiser.h"
#include "types_def.h"

namespace event_model
{
struct TestBatchParams
{
  // dataset folder
  std::string fdir;
  // initial step size of the optimisation
  double iniStep;
  // tolerance that indicates a minimum has been reached
  double tol;
  // maximum iterations
  int maxIter, maxrIter;
  // optimiser status feedback: 0 - no feedback, 1 - feedback at each iteration
  int verbosity;
};

template <typename Dispersion>
int
testBatchExample(const TestBatchParams& testBatchParams)
{
  typedef typename Dispersion::T T;
  typedef typename Dispersion::Model Model;

  constexpr int NDims = Model::NDims, NVars = Model::NVars;

  // optimiser
  typedef GSLfdfOptimiser<Dispersion> Optimiser;

  IO_STATUS ioStatus;
  // read events from file
  Events<T, NDims> evs;
  const std::string fevents(testBatchParams.fdir + "/events.txt");
  ioStatus = load<T, NDims>(fevents, evs);
  if (ioStatus != IO_SUCCESS)
  {
    ioStatusMessage(ioStatus, fevents);
    return -1;
  }
  const int nEvents = evs.size();
  std::cout << "number of events: " << nEvents << '\n';

  int width, height;
  Matrix<T, 3, 3> camParams;
  const std::string fcalib(testBatchParams.fdir + "/calib.txt");
  ioStatus = loadCamParams<T>(fcalib, width, height, camParams);
  if (ioStatus != IO_SUCCESS)
  {
    ioStatusMessage(ioStatus, fcalib);
    return -1;
  }

  Matrix<T> c;
  Vector<T> ts;
  Vector<int> polarity;
  events2eigen<T>(evs, c, ts, polarity);

  // show original events
  showGray(c.template topRows<2>(), polarity, width, height, "Original Events");
  std::cout << "press any key to continue...\n";
  cv::waitKey(0);

  Matrix<T> ct(NDims, nEvents);
  unprojectEvents<T, NDims>()(camParams, c, ct);

  Dispersion dispersion(width);
  dispersion.assignPoints(ct, ts, polarity);

  // initial parameters
  Vector<T, NVars> vars;
  vars.setConstant(1.0e-6);

  // optimise
  Optimiser optimiser(
      dispersion,
      typename Optimiser::OptimiserParams(
          gsl_multimin_fdfminimizer_conjugate_fr, testBatchParams.iniStep,
          testBatchParams.tol, testBatchParams.maxIter,
          testBatchParams.maxrIter, testBatchParams.verbosity));

  std::cout << "score: " << optimiser.run(vars)
            << ", v: " << optimiser.vars().transpose() << '\n';

  // show transformed events according to the estimated parameters
  const Model model;
  Matrix<T> cm(NDims, nEvents), ctm(NDims, nEvents);
  model(optimiser.vars(), ct, ts, ctm);
  projectEvents<T, NDims>()(camParams, ctm, cm);
  showGray(cm.template topRows<2>(), polarity, width, height,
           "Transformed Events");
  std::cout << "press any key to exit...\n";
  cv::waitKey(0);

  return 0;
}
}  // namespace event_model
