#ifndef EVENT_IO_H
#define EVENT_IO_H

#include <Eigen/Core>
#include <fstream>
#include <string>

#include "event/event_type.h"

namespace event
{

template<typename T, int N=Eigen::Dynamic>
bool
load(
  const std::string& fname,
  Events<T, N>& evs)
{
  std::ifstream fin(fname.c_str());
  if(!fin.is_open())
  {
    return false;
  }

  evs.clear();
  Event<T, N> ev;
  while(fin>>ev)
  {
    evs.push_back(ev);
  }

  return true;
}

template<typename T, int N=Eigen::Dynamic>
bool
load(
  const std::string& fname,
  const T& sT,
  const T& eT,
  Events<T, N>& evs)
{
  std::ifstream fin(fname.c_str());
  if(!fin.is_open())
  {
    return false;
  }

  evs.clear();
  Event<T, N> ev;
  while(fin>>ev && ev.ts<sT)
  {}
  if(fin.eof())
  {
    return false;
  }
  do
  {
    evs.push_back(ev);
  } while(fin>>ev && ev.ts<eT);

  return true;
}

template<typename T, int N=Eigen::Dynamic>
bool
load(
  const std::string& fname,
  const T& sT,
  const int nEvents,
  Events<T, N>& evs)
{
  std::ifstream fin(fname.c_str());
  if(!fin.is_open())
  {
    return false;
  }

  evs.clear();
  evs.reserve(nEvents);
  Event<T, N> ev;
  while(fin>>ev && ev.ts<sT)
  {}
  if(fin.eof())
  {
    return false;
  }
  do
  {
    evs.push_back(ev);
  } while(fin>>ev && static_cast<int>(evs.size())<nEvents);

  return true;
}

template<typename T, int N=Eigen::Dynamic>
void
load(
  const int nEvents,
  std::ifstream& fin,
  Events<T, N>& evs)
{
  evs.clear();
  evs.reserve(nEvents);

  Event<T, N> ev;
  while(static_cast<int>(evs.size())<nEvents && fin>>ev)
  {
    evs.push_back(ev);
  }
}

template<typename T, int N>
void
load(
  const int nEvents,
  std::ifstream& fin,
  mtx<T>& c,
  vec<T>& ts,
  vec<int>& polarity)
{
  c.resize(N, nEvents);
  ts.resize(nEvents);
  polarity.resize(nEvents);

  int n;
  Event<T, N> ev;
  for(n=0; n<nEvents && fin>>ev; ++n)
  {
    c.col(n)<<ev.c;
    ts(n)=ev.ts;
    polarity(n)=ev.polarity;
  }

  c.conservativeResize(N, n);
  ts.conservativeResize(n);
  polarity.conservativeResize(n);
}

template<typename T, int N=Eigen::Dynamic>
bool
save(
  const std::string& fname,
  const Events<T, N>& events)
{
  if(events.size()==0)
  {
    return false;
  }

  std::ofstream fout(fname.c_str());
  if(!fout.is_open())
  {
    return false;
  }

  for(int k=0; k<static_cast<int>(events.size()); ++k)
  {
    fout<<events.at(k)<<'\n';
  }

  return true;
}

template<typename T>
bool
save(
  const std::string& fname,
  const mtx<T>& c,
  const vec<T>& ts,
  const vec<int>& polarity)
{
  const int nEvents=c.cols();
  assert(nEvents==ts.size() && nEvents==polarity.size());
  if(nEvents==0)
  {
    return false;
  }

  std::ofstream fout(fname.c_str());
  if(!fout.is_open())
  {
    return false;
  }

  for(int k=0; k<nEvents; ++k)
  {
    fout<<ts(k)<<' '<<c.col(k).transpose()<<' '<<polarity(k)<<'\n';
  }

  return true;
}

}

#endif // EVENT_IO_H
