#include <fstream>
#include <iostream>
#include <string>

#include "EventEMin.h"

using namespace EventEMin;

int
main(int argc, char* argv[])
{
  if (argc < 5)
  {
    std::cout << "usage: " << argv[0]
              << " [events dir] [number of events] "
                 "[saving dir] [file name]\n";
    return -1;
  }

  typedef float T;

  /* you can modify the model used by uncommenting the corresponding line */

  // model
  // typedef IncrementalAffinity<T> Model;
  // typedef IncrementalIsometry<T> Model;
  typedef IncrementalRotation<T> Model;
  // typedef IncrementalSimilarity<T> Model;
  // typedef IncrementalTranslation2D<T> Model;
  // typedef IncrementalTranslationNormal<T> Model;

  constexpr int NDims = Model::NDims;

  /* you can modify the dispersion measure used by uncommenting the
   corresponding line */

  // incremental measures
  typedef IncrementalPotential<Model> Dispersion;
  // typedef IncrementalTsallis<Model> Dispersion;

  // read distorted events from file
  const std::string fevents(std::string(argv[1]) + "/events.txt");
  std::ifstream fin(fevents.c_str());
  if (!fin.is_open())
  {
    std::cerr << "error reading events from file " << fevents << '\n';
    return -1;
  }

  // write estimates to file
  const std::string festimates(std::string(argv[3]) + "/" +
                               std::string(argv[4]) + "_estimates.txt");
  std::ofstream fout(festimates.c_str());
  if (!fout.is_open())
  {
    std::cerr << "error writing estimates to file " << festimates << '\n';
    return -1;
  }

  int width, height;
  cv::Mat camParamsCV, distCoeffs;
  const std::string fcalib(std::string(argv[1]) + "/calib.txt");
  const IO_STATUS ioStatus =
      loadCamParams<T>(fcalib, width, height, camParamsCV, distCoeffs);
  if (ioStatus != IO_SUCCESS)
  {
    ioStatusMessage(ioStatus, fcalib);
    return -1;
  }

  // initialise grid-based undistortion map
  Matrix<T, 3, 3> camParams;
  cv::cv2eigen(camParamsCV, camParams);
  cv::Mat undistortionMap;
  initUndistort<T>(width, height, camParamsCV, distCoeffs, undistortionMap);

  const Vector<T, NDims> scale(Vector<T, NDims>::Ones());

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

  Vector<T, NDims> c, ct;
  T ts;
  int polarity;

  for (int k = 0; undistort<T, NDims>(0, width, 0, height, undistortionMap, fin,
                                      c, ts, polarity) == IO_SUCCESS;
       ++k)
  {
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
