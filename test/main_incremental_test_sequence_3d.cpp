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
                 "[min depth] [depth scale] [saving dir] [file "
                 "name]\n";
    return -1;
  }

  typedef float T;

  /* you can modify the model used by uncommenting the corresponding line */

  // model
  typedef IncrementalSixDOF<T> Model;
  // typedef IncrementalTranslation3D<T> Model;

  constexpr int NDims = Model::NDims;

  /* you can modify the dispersion measure used by uncommenting the
   corresponding line */

  // incremental measures
  typedef IncrementalPotential<Model> Dispersion;
  // typedef IncrementalTsallis<Model> Dispersion;

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

  Vector<T, NDims> scale;
  scale << T(1.0), T(1.0), std::atof(argv[4]);

  // tolerance that indicates a minimum has been reached
  const T minStep = T(1.0e-6);
  // maximum iterations
  const int maxIter = 10;
  // neighbouring radius
  const int wSize = 4;
  // number of events to maintain
  const int nEvents = std::atoi(argv[2]);
  Dispersion dispersion(camParams, scale,
                        Dispersion::Params(minStep, maxIter, wSize), nEvents,
                        {width, height});

  const T depthThresh = std::atof(argv[3]);
  Vector<T, NDims> c, ct;
  T ts;
  int polarity;

  for (int k = 0; load<T, NDims>(fin, c, ts, polarity) == IO_SUCCESS; ++k)
  {
    const int x = std::round(c(0));
    const int y = std::round(c(1));

    if (!(0 <= x && x < width && 0 <= y && y < height) || c(2) < depthThresh)
    {
      --k;
      continue;
    }
    unprojectEvent<T, NDims>()(camParams, c, ct);

    dispersion.run(ct, ts);

    if ((k + 1) % nEvents == 0)
    {
      std::cout << "ts: " << ts << ", vars: " << dispersion.vars().transpose()
                << '\n';
      fout << ts << ' ' << dispersion.vars().transpose() << std::endl;
    }
  }

  return 0;
}
