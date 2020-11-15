#include <fstream>
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
main(int argc, char* argv[])
{
  if(argc<5)
  {
    std::cout<<"usage: "<<argv[0]<<" [events dir] [number of events] [saving dir] [file name]\n";
    return -1;
  }

  typedef float T;

  // model
  typedef model::Rotation<T> Model;

  constexpr int ndims=Model::ndims(), nvars=Model::nvars();

  // you can modify the entropy used by uncommenting the corresponding line
  //typedef dispersion::ApproxPotential<Model> Dispersion;
  //typedef dispersion::ApproxRenyi<Model> Dispersion;
  //typedef dispersion::ApproxShannon<Model> Dispersion;
  //typedef dispersion::ApproxSharmaMittal<Model> Dispersion;
  typedef dispersion::ApproxTsallis<Model> Dispersion;

  // optimiser
  typedef optimiser::GSLfdfOptimiser<Dispersion> Optimiser;

  // read distorted events from file
  const std::string fevents(std::string(argv[1])+"/events.txt");
  std::ifstream fin(fevents.c_str());
  if(!fin.is_open())
  {
    std::cerr<<"error reading events from file "<<fevents<<'\n';
    return -1;
  }

  // write estimates to file
  const std::string festimates(std::string(argv[3])+"/"+std::string(argv[4])+"_estimates.txt");
  std::ofstream fout(festimates.c_str());
  if(!fout.is_open())
  {
    std::cerr<<"error writing estimates to file "<<festimates<<'\n';
    return -1;
  }

  int width, height;
  cv::Mat camParamsCV, distCoeffs;
  const std::string fcalib(std::string(argv[1])+"/calib.txt");
  if(!loadCamParams<T>(fcalib, width, height, camParamsCV, distCoeffs))
  {
    std::cout<<"error reading calibration from file "<<fcalib<<'\n';
    return -1;
  }

  // initialise grid-based undistortion map
  mtx<T, 3, 3> camParams;
  cv::cv2eigen(camParamsCV, camParams);
  cv::Mat undistortionMap;
  initUndistort<T>(width, height, camParamsCV, distCoeffs, undistortionMap);

  Dispersion dispersion({width, height});

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

  const int nEvents=std::atoi(argv[2]);
  mtx<T> c;
  vec<T> ts;
  vec<int> polarity;
  while(!fin.eof())
  {
    // load batch of events and undistort them
    event::undistort<T>(0, width-1, 0, height-1,
      undistortionMap, nEvents, fin, c, ts, polarity);
    mtx<T> ct(ndims, c.cols());
    event::transformEvents<T>(camParams, c, ct);
    dispersion.assignPoints(ct, ts, polarity);

    // optimise
    Optimiser optimiser(dispersion, Optimiser::OptimiserParams(
      gsl_multimin_fdfminimizer_conjugate_fr, iniStep, tol, maxIter, maxrIter, verbosity));

    optimiser.run(vars);
    vars=optimiser.vars();
    std::cout<<"ts: "<<dispersion.tsend()<<", vars: "<<vars.transpose()<<'\n';
    fout<<dispersion.tsend()<<' '<<vars.transpose()<<std::endl;
  }

  return 0;
}
