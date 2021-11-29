#include "EventEMin.h"

using namespace EventEMin;

int
main(void)
{
  typedef float T;

  // model
  typedef IncrementalSixDOF<T> Model;

  /* you can modify the dispersion measure used by uncommenting the
   corresponding line */

  // incremental measures
  typedef IncrementalPotential<Model> Dispersion;
  // typedef IncrementalTsallis<Model> Dispersion;

  // test parameters for incremental mode - default parameters
  TestIncrementalParams testIncrementalParams;
  // dataset folder
  testIncrementalParams.fdir = "../dataset/indoor_flying1";
  // tolerance that indicates a minimum has been reached
  testIncrementalParams.minStep = T(1.0e-6);
  // maximum iterations
  testIncrementalParams.maxIter = 10;
  // neighbouring radius
  testIncrementalParams.wSize = 4;
  // number of events to maintain
  testIncrementalParams.nEvents = 15000;

  return testIncrementalExample<Dispersion>(testIncrementalParams);
}
