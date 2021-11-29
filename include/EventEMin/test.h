#ifndef EVENT_EMIN_TEST_H
#define EVENT_EMIN_TEST_H

#include <iostream>
#include <string>

#include "EventEMin/event.h"
#include "EventEMin/image.h"
#include "EventEMin/optimiser.h"
#include "EventEMin/types_def.h"

namespace EventEMin
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

struct TestIncrementalParams
{
  // dataset folder
  std::string fdir;
  // tolerance that indicates a minimum has been reached
  double minStep;
  // maximum iterations
  int maxIter;
  // neighbouring radius
  int wSize;
  // number of events to maintain
  int nEvents;
};

template <typename Dispersion>
int
testIncrementalExample(const TestIncrementalParams& testIncrementalParams)
{
  typedef typename Dispersion::T T;
  typedef typename Dispersion::Model Model;

  constexpr int NDims = Model::NDims;

  IO_STATUS ioStatus;
  // read events from file
  Events<T, NDims> evs;
  const std::string fevents(testIncrementalParams.fdir + "/events.txt");
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
  const std::string fcalib(testIncrementalParams.fdir + "/calib.txt");
  ioStatus = loadCamParams<T>(fcalib, width, height, camParams);
  if (ioStatus != IO_SUCCESS)
  {
    ioStatusMessage(ioStatus, fcalib);
    return -1;
  }

  const std::string imgTransformedEvents("Transformed Events"),
      imgOriginalEvents("Original Events");
  cv::namedWindow(imgTransformedEvents);
  cv::moveWindow(imgTransformedEvents, 440, 0);
  cv::namedWindow(imgOriginalEvents);

  Matrix<T> c;
  Vector<T> ts;
  Vector<int> polarity;
  events2eigen<T>(evs, c, ts, polarity);

  Matrix<T> ct(NDims, nEvents);
  unprojectEvents<T, NDims>()(camParams, c, ct);

  Vector<T, NDims> scale;
  scale.template head<NDims>() = camParams.diagonal().template head<NDims>();
  if (NDims == 3)
  {
    scale(2) = T(100.0);
  }

  Dispersion dispersion(
      camParams, scale,
      typename Dispersion::Params(testIncrementalParams.minStep,
                                  testIncrementalParams.maxIter,
                                  testIncrementalParams.wSize),
      testIncrementalParams.nEvents, {width, height});

  for (int k = 0; k < nEvents; ++k)
  {
    const int x = std::round(c(0, k));
    const int y = std::round(c(1, k));
    if (0 <= x && x < width && 0 <= y && y < height)
    {
      dispersion.run(ct.col(k), ts(k));
    }

    if ((k + 1) % testIncrementalParams.nEvents == 0)
    {
      std::cout << "ts: " << ts(k) << ", v: " << -dispersion.vars().transpose()
                << '\n';

      const int kk = k + 1 - testIncrementalParams.nEvents;
      const Model model;
      Matrix<T> cm(NDims, testIncrementalParams.nEvents),
          ctm(NDims, testIncrementalParams.nEvents);
      model(dispersion.vars(), ct.middleCols(kk, testIncrementalParams.nEvents),
            ts.segment(kk, testIncrementalParams.nEvents), ctm);
      projectEvents<T, NDims>()(camParams, ctm, cm);

      showGray(cm.template topRows<2>(),
               polarity.segment(kk, testIncrementalParams.nEvents), width,
               height, imgTransformedEvents);
      showGray(
          c.middleCols(kk, testIncrementalParams.nEvents).template topRows<2>(),
          polarity.segment(kk, testIncrementalParams.nEvents), width, height,
          imgOriginalEvents);
      cv::waitKey(30);
    }
  }

  std::cout << "v: " << -dispersion.vars().transpose() << '\n';
  std::cout << "press any key to exit...\n";
  cv::waitKey(0);

  return 0;
}
}  // namespace EventEMin

#endif  // EVENT_EMIN_TEST_H
