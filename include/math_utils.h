#ifndef MATH_UTILS
#define MATH_UTILS

#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <utility>
#include <algorithm> // for std::minmax_element

// #include <math.h>
#define M_PI 3.14159265358979323846264338327950288 /* pi             */

namespace Geo {

// Hash function
struct hashPairIntFunction {
  size_t operator()(const std::pair<int, int> &x) const {
    int h1 = std::hash<int>{}(x.first);
    int h2 = std::hash<int>{}(x.second);
    return h1 ^ h2;
  }
};

template <class ScalarType> inline ScalarType ToDeg(const ScalarType &a) {
  return a * (ScalarType)180.0f / ScalarType(M_PI);
}

template <class ScalarType> inline ScalarType normalizeAngle(ScalarType angle) {
  angle = (ScalarType)fmod((float)angle, (float)2 * M_PI);
  if (angle < 0) {
    angle += 2 * M_PI;
  }
  return angle;
}

template <class ScalarType>
inline ScalarType MaxPercentile(const std::vector<ScalarType> &Values,
                                const ScalarType &Percentile) {
  assert((Percentile >= 0) && (Percentile <= 1));

  // then compute the max (considering the percentile)
  std::vector<ScalarType> ValuesSorted = Values;
  std::sort(ValuesSorted.begin(), ValuesSorted.end());
  int IndexMax = floor(
      ((ScalarType)ValuesSorted.size() * ((ScalarType)1 - Percentile)) + 0.5);
  IndexMax = std::min(IndexMax, (int)ValuesSorted.size() - 1);
  IndexMax = std::max(IndexMax, 0);
  ScalarType TestVal = ValuesSorted[IndexMax];
  return TestVal;
}

template <class ScalarType>
inline void Trim(std::vector<ScalarType> &Values, const ScalarType &MaxVal) {
  for (size_t i = 0; i < Values.size(); i++)
    Values[i] = std::min(Values[i], MaxVal);
}

template <class ScalarType>
inline ScalarType LinearMap(ScalarType &Value, 
                            const ScalarType Min0,
                            const ScalarType Max0, 
                            const ScalarType Min1,
                            const ScalarType Max1) {
  assert(Max0 > Min0);
  Value = std::max(Value, Min0);
  Value = std::min(Value, Max0);

  assert(Max1 > Min1);
  ScalarType ratio = (Value - Min0) / (Max0 - Min0);
  return (Min1 + (Max1 - Min1) * ratio);
}

template <class ScalarType>
inline void LinearMap(std::vector<ScalarType> &Values, 
                     const ScalarType Min1,
                     const ScalarType Max1) 
{
  assert(Max1 > Min1);
  
  ScalarType MinV0=std::numeric_limits<ScalarType>::max();
  ScalarType MaxV0=std::numeric_limits<ScalarType>::lowest();
  for (size_t i = 0; i < Values.size(); i++)
  {
    MinV0=std::min(MinV0,Values[i]);
    MaxV0=std::max(MaxV0,Values[i]);
  }

  for (size_t i = 0; i < Values.size(); i++)
    Values[i] = LinearMap(Values[i], MinV0,MaxV0, Min1, Max1);
}

template <class ScalarType>
inline void TrimPercentile(std::vector<ScalarType> &Values, 
                     const ScalarType minP,
                     const ScalarType maxP) 
{
  assert(minP >=0);
  assert(minP <=1);
  assert(maxP >=0);
  assert(maxP <=1);

  std::vector<ScalarType> ValuesSorted = Values;
  std::sort(ValuesSorted.begin(), ValuesSorted.end());
  int IndexMax = floor(((ScalarType)ValuesSorted.size() * ((ScalarType)1 - maxP)) + 0.5);
  int IndexMin = floor(((ScalarType)ValuesSorted.size() * ((ScalarType)minP)) + 0.5);

  IndexMax = std::min(IndexMax, (int)ValuesSorted.size() - 1);
  IndexMax = std::max(IndexMax, 0);
  IndexMin = std::min(IndexMin, (int)ValuesSorted.size() - 1);
  IndexMin = std::max(IndexMin, 0);

  ScalarType ValueMin=ValuesSorted[IndexMin];
  ScalarType ValueMax=ValuesSorted[IndexMax];

  for (size_t i = 0; i < Values.size(); i++)
  {
    Values[i] = std::max(Values[i],ValueMin);
    Values[i] = std::min(Values[i],ValueMax);
  }
}
}

#endif
