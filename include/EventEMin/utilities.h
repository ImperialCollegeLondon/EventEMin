#ifndef EVENT_EMIN_UTILITIES_H
#define EVENT_EMIN_UTILITIES_H

#include <cmath>

namespace EventEMin
{
inline double
computeExp(const double x)
{
#ifdef EventEMin_FAST_EXP
  // fast computation of exp(x)
  if (x < -100)
  {
    return 0;
  }
  union
  {
    double d;
    long l;
  } u;
  u.l = long(double((long(1) << 52) / 0.69314718056) * x +
             double((long(1) << 52) * 1023));
  return u.d;
#else
  return std::exp(x);
#endif
}

// add modulus (cyclic)
inline int
addCyclic(const int x, const int a, const int n)
{
  return (x + a + n) % n;
}
// decrement modulus (cyclic)
inline int
decrementCyclic(const int x, const int n)
{
  return addCyclic(x, -1, n);
}
// increment modulus (cyclic)
inline int
incrementCyclic(const int x, const int n)
{
  return addCyclic(x, 1, n);
}

// fast cyclic operations assume n = 2^m-1, where m > 0 is some integer
// fast add modulus (cyclic)
inline int
fastAddCyclic(const int x, const int a, const int n)
{
  assert(0 < n && ((n + 1) & n) == 0);
  return (x + a) & n;
}
// fast decrement modulus (cyclic)
inline int
fastDecrementCyclic(const int x, const int n)
{
  return fastAddCyclic(x, -1, n);
}

// fast increment modulus (cyclic)
inline int
fastIncrementCyclic(const int x, const int n)
{
  return fastAddCyclic(x, 1, n);
}

// computes the next power of 2 higher than n
inline int
nextPower2(int n)
{
  int i = 0;
  for (--n; n > 0; n >>= 1, ++i)
    ;
  return 1 << i;
}
}  // namespace EventEMin

#endif  // EVENT_EMIN_UTILITIES_H
