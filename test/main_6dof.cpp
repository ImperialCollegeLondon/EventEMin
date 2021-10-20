#include "dispersion.h"
#include "model.h"
#include "test.h"

using namespace event_model;

int
main(void)
{
  typedef float T;

  // model
  typedef SixDOF<T> Model;

  /* you can modify the dispersion measure used by uncommenting the
   corresponding line */

  // exact measures
  // typedef Potential<Model> Dispersion;
  // typedef Renyi<Model> Dispersion;
  // typedef Shannon<Model> Dispersion;
  // typedef SharmaMittal<Model> Dispersion;
  // typedef Tsallis<Model> Dispersion;
  // approximate measures
  // typedef ApproxPotential<Model> Dispersion;
  // typedef ApproxRenyi<Model> Dispersion;
  // typedef ApproxShannon<Model> Dispersion;
  // typedef ApproxSharmaMittal<Model> Dispersion;
  typedef ApproxTsallis<Model> Dispersion;

  // test parameters for batch mode - default parameters
  TestBatchParams testBatchParams;
  // dataset folder
  testBatchParams.fdir = "../dataset/indoor_flying1";
  // initial step size of the optimisation
  testBatchParams.iniStep = 1.0;
  // tolerance that indicates a minimum has been reached
  testBatchParams.tol = 1.0e-16;
  // maximum iterations
  testBatchParams.maxIter = 100, testBatchParams.maxrIter = 1;
  // optimiser status feedback: 0 - no feedback, 1 - feedback at each iteration
  testBatchParams.verbosity = 1;

  return testBatchExample<Dispersion>(testBatchParams);
}
