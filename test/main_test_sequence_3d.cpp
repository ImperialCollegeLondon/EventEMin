#include <fstream>
#include <iostream>
#include <string>

#include "EventEMin.h"

using namespace EventEMin;

int
main(int argc, char* argv[])
{
  if (argc < 7)
  {
    std::cout << "usage: " << argv[0]
              << " [events dir] [number of events] "
                 "[min depth] [max depth] [saving dir] "
                 "[file name]\n";
    return -1;
  }

  typedef float T;

  /* you can modify the model used by uncommenting the corresponding line */

  // model
  typedef SixDOF<T> Model;
  // typedef Translation3D<T> Model;

  constexpr int NDims = Model::NDims, NVars = Model::NVars;

  // exact measures
  // typedef Potential<Model> Dispersion;
  // typedef Renyi<Model> Dispersion;
  // typedef Shannon<Model> Dispersion;
  // typedef SharmaMittal<Model> Dispersion;
  // typedef Tsallis<Model> Dispersion;
  // approximate measures
  // typedef ApproximatePotential<Model> Dispersion;
  // typedef ApproximateRenyi<Model> Dispersion;
  // typedef ApproximateShannon<Model> Dispersion;
  // typedef ApproximateSharmaMittal<Model> Dispersion;
  typedef ApproximateTsallis<Model> Dispersion;

  // optimiser
  typedef GSLfdfOptimiser<Dispersion> Optimiser;

  // read undistorted depth-augmented events from file
  const std::string fevents(std::string(argv[1]) + "/events.txt");
  std::ifstream fin(fevents.c_str());
  if (!fin.is_open())
  {
    std::cerr << "error reading events from file " << fevents << '\n';
    return -1;
  }

  // write estimates to file
  const std::string festimates(std::string(argv[5]) + "/" +
                               std::string(argv[6]) + "_estimates.txt");
  std::ofstream fout(festimates.c_str());
  if (!fout.is_open())
  {
    std::cerr << "error writing estimates to file " << festimates << '\n';
    return -1;
  }

  int width, height;
  Matrix<T, 3, 3> camParams;
  const std::string fcalib(std::string(argv[1]) + "/calib.txt");
  const IO_STATUS ioStatus = loadCamParams<T>(fcalib, width, height, camParams);
  if (ioStatus != IO_SUCCESS)
  {
    ioStatusMessage(ioStatus, fcalib);
    return -1;
  }

  Dispersion dispersion(width);

  // initial parameters
  Vector<T, NVars> vars;
  vars.setConstant(1.0e-6);

  // initial step size of the optimisation
  const double iniStep = 1.0;
  // tolerance that indicates a minimum has been reached
  const double tol = 1.0e-16;
  // maximum iterations
  const int maxIter = 100, maxrIter = 1;
  // optimiser status feedback: 0 - no feedback, 1 - feedback at each iteration
  const int verbosity = 1;
  // apply whitening pre-processing step
  const bool whiten = false;

  const int nEvents = std::atoi(argv[2]);
  Matrix<T> c;
  Vector<T> ts;
  Vector<int> polarity;
  while (loadDepthThresh<T>(nEvents, std::atof(argv[3]), std::atof(argv[4]),
                            fin, c, ts, polarity) == IO_SUCCESS)
  {
    Matrix<T> ct(NDims, c.cols());
    unprojectEvents<T, NDims>()(camParams, c, ct);
    dispersion.assignPoints(ct, ts, polarity, whiten);

    Optimiser optimiser(
        dispersion,
        Optimiser::OptimiserParams(gsl_multimin_fdfminimizer_conjugate_fr,
                                   iniStep, tol, maxIter, maxrIter, verbosity));

    optimiser.run(vars);
    vars = optimiser.vars();
    std::cout << "ts: " << dispersion.tsEnd() << ", vars: " << vars.transpose()
              << '\n';
    fout << dispersion.tsEnd() << ' ' << vars.transpose() << std::endl;
  }

  return 0;
}
