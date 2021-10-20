#ifndef EVENT_IO_H
#define EVENT_IO_H

#include <fstream>
#include <iostream>
#include <string>

#include "event/type.h"
#include "types_def.h"

namespace event_model
{
enum IO_STATUS
{
  IO_FAIL = -1,
  IO_SUCCESS,
  IO_EMPTY,
  IO_LESS
};

void
ioStatusMessage(const IO_STATUS status, const std::string& fname = "")
{
  switch (status)
  {
    case IO_FAIL:
      std::cout << "fail opening file " << fname << '\n';
      return;
    case IO_SUCCESS:
      std::cout << "success reading file " << fname << '\n';
      return;
    case IO_EMPTY:
      std::cout << "no events read from file " << fname << '\n';
      return;
    case IO_LESS:
      std::cout << "less events read from file " << fname << '\n';
      return;
  }
}

template <typename T>
IO_STATUS
loadCamParams(const std::string& fname, int& width, int& height,
              CvMatrix& camParams, CvMatrix& distCoeffs)
{
  std::ifstream fin(fname.c_str());
  if (!fin.is_open())
  {
    return IO_FAIL;
  }

  camParams.create(3, 3, CV_TYPE(T, 1));
  distCoeffs.create(5, 1, CV_TYPE(T, 1));
  fin >> camParams.at<T>(0, 0) >> camParams.at<T>(1, 1) >>
      camParams.at<T>(0, 2) >> camParams.at<T>(1, 2);
  camParams.at<T>(1, 1) = -camParams.at<T>(1, 1);
  camParams.at<T>(0, 1) = 0.0;
  camParams.at<T>(1, 0) = 0.0;
  camParams.at<T>(2, 0) = 0.0;
  camParams.at<T>(2, 1) = 0.0;
  camParams.at<T>(2, 2) = 1.0;
  fin >> distCoeffs.at<T>(0) >> distCoeffs.at<T>(1) >> distCoeffs.at<T>(2) >>
      distCoeffs.at<T>(3) >> distCoeffs.at<T>(4);
  fin >> width >> height;

  return IO_SUCCESS;
}

template <typename T>
IO_STATUS
loadCamParams(const std::string& fname, int& width, int& height,
              CvMatrix& camParams)
{
  CvMatrix distCoeffs;
  return loadCamParams<T>(fname, width, height, camParams, distCoeffs);
}

template <typename T>
IO_STATUS
loadCamParams(const std::string& fname, int& width, int& height,
              Matrix<T, 3, 3>& camParams)
{
  CvMatrix camParamsCV;
  const IO_STATUS status = loadCamParams<T>(fname, width, height, camParamsCV);
  if (status == IO_SUCCESS)
  {
    cv::cv2eigen(camParamsCV, camParams);
  }
  return status;
}

template <typename T, int N>
IO_STATUS
load(std::ifstream& fin, Event<T, N>& ev)
{
  if (fin.eof())
  {
    return IO_EMPTY;
  }
  fin >> ev;
  return IO_SUCCESS;
}

template <typename T, int N>
IO_STATUS
load(std::ifstream& fin, Ref<Vector<T, N> > c, T& ts, int& polarity)
{
  Event<T, N> ev;
  const IO_STATUS ioStatus = load<T, N>(fin, ev);
  if (ioStatus == IO_SUCCESS)
  {
    event2eigen(ev, c, ts, polarity);
    return IO_SUCCESS;
  }
  return ioStatus;
}

template <typename T, int N>
IO_STATUS
load(const std::string& fname, Events<T, N>& evs)
{
  std::ifstream fin(fname.c_str());
  if (!fin.is_open())
  {
    return IO_FAIL;
  }

  evs.clear();
  Event<T, N> ev;
  while (load<T, N>(fin, ev) == IO_SUCCESS)
  {
    evs.push_back(ev);
  }

  if (evs.size() > 0)
  {
    return IO_SUCCESS;
  }
  return IO_EMPTY;
}

template <typename T, int N>
IO_STATUS
load(const std::string& fname, const T& sT, const T& eT, Events<T, N>& evs)
{
  std::ifstream fin(fname.c_str());
  if (!fin.is_open())
  {
    return IO_FAIL;
  }

  evs.clear();
  Event<T, N> ev;
  while (load<T, N>(fin, ev) == IO_SUCCESS && ev.ts < sT)
  {
  }
  if (fin.eof())
  {
    return IO_EMPTY;
  }
  do
  {
    evs.push_back(ev);
  } while (load<T, N>(fin, ev) == IO_SUCCESS && ev.ts < eT);

  return IO_SUCCESS;
}

template <typename T, int N>
IO_STATUS
load(const std::string& fname, const T& sT, const int nEvents,
     Events<T, N>& evs)
{
  std::ifstream fin(fname.c_str());
  if (!fin.is_open())
  {
    return IO_FAIL;
  }

  evs.clear();
  evs.reserve(nEvents);
  Event<T, N> ev;
  while (load<T, N>(fin, ev) == IO_SUCCESS && ev.ts < sT)
  {
  }
  if (fin.eof())
  {
    return IO_EMPTY;
  }
  do
  {
    evs.push_back(ev);
  } while (load<T, N>(fin, ev) == IO_SUCCESS &&
           static_cast<int>(evs.size()) < nEvents);

  if (static_cast<int>(evs.size()) == nEvents)
  {
    return IO_SUCCESS;
  }
  return IO_LESS;
}

template <typename T, int N>
IO_STATUS
load(const int nEvents, std::ifstream& fin, Events<T, N>& evs)
{
  if (fin.eof())
  {
    return IO_EMPTY;
  }

  evs.clear();
  evs.reserve(nEvents);

  Event<T, N> ev;
  while (load<T, N>(fin, ev) == IO_SUCCESS &&
         static_cast<int>(evs.size()) < nEvents)
  {
    evs.push_back(ev);
  }

  if (static_cast<int>(evs.size()) == nEvents)
  {
    return IO_SUCCESS;
  }
  return IO_LESS;
}

template <typename T, int N>
IO_STATUS
load(const int nEvents, std::ifstream& fin, Matrix<T>& c, Vector<T>& ts,
     Vector<int>& polarity)
{
  if (fin.eof())
  {
    return IO_EMPTY;
  }

  c.resize(N, nEvents);
  ts.resize(nEvents);
  polarity.resize(nEvents);

  IO_STATUS ioStatus = IO_SUCCESS;
  int n;
  for (n = 0; ioStatus == IO_SUCCESS && n < nEvents; ++n)
  {
    ioStatus = load<T, N>(fin, c.col(n), ts(n), polarity(n));
  }

  if (n == nEvents)
  {
    return IO_SUCCESS;
  }
  c.conservativeResize(N, n);
  ts.conservativeResize(n);
  polarity.conservativeResize(n);
  return IO_LESS;
}

template <typename T, int N>
IO_STATUS
load(const T& eT, std::ifstream& fin, Events<T, N>& evs)
{
  if (fin.eof())
  {
    return IO_EMPTY;
  }

  evs.clear();
  Event<T, N> ev;
  while (load<T, N>(fin, ev) == IO_SUCCESS && ev.ts < eT)
  {
    evs.push_back(ev);
  }

  return IO_SUCCESS;
}

template <typename T, int N>
IO_STATUS
load(const T& eT, std::ifstream& fin, Matrix<T>& c, Vector<T>& ts,
     Vector<int>& polarity)
{
  Events<T, N> evs;
  if (load<T, N>(eT, fin, evs) == IO_EMPTY)
  {
    return IO_EMPTY;
  }

  events2eigen(evs, c, ts, polarity);
  return IO_SUCCESS;
}

template <typename T>
IO_STATUS
loadDepthThresh(const std::string& fname, const T& sT, const int nEvents,
                const T& depthMin, const T& depthMax, Events<T, 3>& evs)
{
  std::ifstream fin(fname.c_str());
  if (!fin.is_open())
  {
    return IO_FAIL;
  }

  evs.clear();
  evs.reserve(nEvents);
  Event<T, 3> ev;

  while (load<T, 3>(fin, ev) == IO_SUCCESS && ev.ts < sT)
  {
  }
  if (fin.eof())
  {
    return IO_EMPTY;
  }
  do
  {
    if (!(ev.c(2) < depthMin || ev.c(2) > depthMax))
    {
      evs.push_back(ev);
    }
  } while (load<T, 3>(fin, ev) == IO_SUCCESS &&
           static_cast<int>(evs.size()) < nEvents);

  if (static_cast<int>(evs.size()) == nEvents)
  {
    return IO_SUCCESS;
  }
  return IO_LESS;
}

template <typename T>
IO_STATUS
loadDepthThresh(const int nEvents, const T& depthMin, const T& depthMax,
                std::ifstream& fin, Matrix<T>& c, Vector<T>& ts,
                Vector<int>& polarity)
{
  if (fin.eof())
  {
    return IO_EMPTY;
  }

  c.resize(3, nEvents);
  ts.resize(nEvents);
  polarity.resize(nEvents);

  IO_STATUS ioStatus = IO_SUCCESS;
  int n;
  for (n = 0; ioStatus == IO_SUCCESS && n < nEvents; ++n)
  {
    ioStatus = load<T, 3>(fin, c.col(n), ts(n), polarity(n));
    if (ioStatus == IO_SUCCESS && (c(2, n) < depthMin || c(2, n) > depthMax))
    {
      --n;
    }
  }

  if (n == nEvents)
  {
    return IO_SUCCESS;
  }
  c.conservativeResize(3, n);
  ts.conservativeResize(n);
  polarity.conservativeResize(n);
  return IO_LESS;
}

template <typename T, int N>
IO_STATUS
save(const std::string& fname, const Events<T, N>& evs)
{
  if (evs.size() == 0)
  {
    return IO_EMPTY;
  }

  std::ofstream fout(fname.c_str());
  if (!fout.is_open())
  {
    return IO_FAIL;
  }

  for (int k = 0; k < static_cast<int>(evs.size()); ++k)
  {
    fout << evs[k] << '\n';
  }

  return IO_SUCCESS;
}

template <typename T>
IO_STATUS
save(const std::string& fname, const Matrix<T>& c, const Vector<T>& ts,
     const Vector<int>& polarity)
{
  const int nEvents = c.cols();
  assert(ts.size() == nEvents);
  assert(polarity.size() == nEvents);

  if (nEvents == 0)
  {
    return IO_EMPTY;
  }

  std::ofstream fout(fname.c_str());
  if (!fout.is_open())
  {
    return IO_FAIL;
  }

  for (int k = 0; k < nEvents; ++k)
  {
    fout << ts(k) << ' ' << c.col(k).transpose() << ' ' << polarity(k) << '\n';
  }

  return IO_SUCCESS;
}

template <typename T, int N>
IO_STATUS
undistort(const int xMin, const int xMax, const int yMin, const int yMax,
          const CvMatrix& map, std::ifstream& fin, Event<T, N>& ev)
{
  assert(N >= 2);
  while (!fin.eof())
  {
    fin >> ev;
    CvVector<T, 2> mapPoint =
        map.at<CvVector<T, 2> >(static_cast<int>(std::round(ev.c(1))),
                                static_cast<int>(std::round(ev.c(0))));
    const int x = static_cast<int>(std::round(mapPoint[0]));
    const int y = static_cast<int>(std::round(mapPoint[1]));
    if (xMin <= x && x < xMax && yMin <= y && y < yMax)
    {
      ev.c.template head<2>() << mapPoint[0], mapPoint[1];
      return IO_SUCCESS;
    }
  }
  return IO_EMPTY;
}

template <typename T, int N>
IO_STATUS
undistort(const int xMin, const int xMax, const int yMin, const int yMax,
          const CvMatrix& map, std::ifstream& fin, Event<T, N>& ev,
          Event<T, N>& evo)
{
  assert(N >= 2);
  while (!fin.eof())
  {
    fin >> evo;
    CvVector<T, 2> mapPoint =
        map.at<CvVector<T, 2> >(static_cast<int>(std::round(evo.c(1))),
                                static_cast<int>(std::round(evo.c(0))));
    const int x = static_cast<int>(std::round(mapPoint[0]));
    const int y = static_cast<int>(std::round(mapPoint[1]));
    if (xMin <= x && x < xMax && yMin <= y && y < yMax)
    {
      ev.c.template head<2>() << mapPoint[0], mapPoint[1];
      ev.ts = evo.ts;
      ev.polarity = evo.polarity;
      return IO_SUCCESS;
    }
  }
  return IO_EMPTY;
}

template <typename T, int N>
IO_STATUS
undistort(const int xMin, const int xMax, const int yMin, const int yMax,
          const CvMatrix& map, std::ifstream& fin, Ref<Vector<T, N> > c, T& ts,
          int& polarity)
{
  Event<T, N> ev;
  const IO_STATUS ioStatus =
      undistort<T, N>(xMin, xMax, yMin, yMax, map, fin, ev);
  if (ioStatus == IO_SUCCESS)
  {
    event2eigen(ev, c, ts, polarity);
  }
  return ioStatus;
}

template <typename T, int N>
IO_STATUS
undistort(const int xMin, const int xMax, const int yMin, const int yMax,
          const CvMatrix& map, std::ifstream& fin, Ref<Vector<T, N> > c, T& ts,
          int& polarity, Ref<Vector<T, N> > co)
{
  Event<T, N> ev, evo;
  const IO_STATUS ioStatus =
      undistort<T, N>(xMin, xMax, yMin, yMax, map, fin, ev, evo);
  if (ioStatus == IO_SUCCESS)
  {
    event2eigen(ev, c, ts, polarity);
    event2eigen(evo, co, ts, polarity);
  }
  return ioStatus;
}

template <typename T, int N>
IO_STATUS
undistort(const int xMin, const int xMax, const int yMin, const int yMax,
          const CvMatrix& map, const int nEvents, std::ifstream& fin,
          Events<T, N>& evs)
{
  assert(N >= 2);
  if (fin.eof())
  {
    return IO_EMPTY;
  }

  evs.clear();
  evs.reserve(nEvents);

  IO_STATUS ioStatus = IO_SUCCESS;
  Event<T, N> ev;
  while (ioStatus != IO_EMPTY && static_cast<int>(evs.size()) < nEvents)
  {
    ioStatus = undistort<T, N>(xMin, xMax, yMin, yMax, map, fin, ev);
    if (ioStatus == IO_SUCCESS)
    {
      evs.push_back(ev);
    }
  }
  if (static_cast<int>(evs.size()) == nEvents)
  {
    return IO_SUCCESS;
  }
  return IO_LESS;
}

template <typename T, int N>
IO_STATUS
undistort(const int xMin, const int xMax, const int yMin, const int yMax,
          const CvMatrix& map, const int nEvents, std::ifstream& fin,
          Matrix<T>& c, Vector<T>& ts, Vector<int>& polarity)
{
  assert(N >= 2);
  if (fin.eof())
  {
    return IO_EMPTY;
  }

  c.resize(N, nEvents);
  ts.resize(nEvents);
  polarity.resize(nEvents);

  IO_STATUS ioStatus = IO_SUCCESS;
  int n = 0;
  while (ioStatus != IO_EMPTY && n < nEvents)
  {
    ioStatus = undistort<T, N>(xMin, xMax, yMin, yMax, map, fin, c.col(n),
                               ts(n), polarity(n));
    if (ioStatus == IO_SUCCESS)
    {
      ++n;
    }
  }

  if (n == nEvents)
  {
    return IO_SUCCESS;
  }
  c.conservativeResize(N, n);
  ts.conservativeResize(n);
  polarity.conservativeResize(n);
  return IO_LESS;
}

template <typename T, int N>
IO_STATUS
undistort(const int xMin, const int xMax, const int yMin, const int yMax,
          const CvMatrix& map, const int nEvents, std::ifstream& fin,
          Matrix<T>& c, Vector<T>& ts, Vector<int>& polarity, Matrix<T>& co)
{
  assert(N >= 2);
  if (fin.eof())
  {
    return IO_EMPTY;
  }

  c.resize(N, nEvents);
  ts.resize(nEvents);
  polarity.resize(nEvents);
  co.resize(N, nEvents);

  IO_STATUS ioStatus = IO_SUCCESS;
  int n = 0;
  while (ioStatus != IO_EMPTY && n < nEvents)
  {
    ioStatus = undistort<T, N>(xMin, xMax, yMin, yMax, map, fin, c.col(n),
                               ts(n), polarity(n), co.col(n));
    if (ioStatus == IO_SUCCESS)
    {
      ++n;
    }
  }

  if (n == nEvents)
  {
    return IO_SUCCESS;
  }
  c.conservativeResize(N, n);
  ts.conservativeResize(n);
  polarity.conservativeResize(n);
  co.conservativeResize(N, n);
  return IO_LESS;
}
}  // namespace event_model

#endif  // EVENT_IO_H
