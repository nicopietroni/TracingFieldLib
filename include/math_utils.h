#ifndef MATH_UTILS
#define MATH_UTILS

#include <cmath>
#include <utility>
#include <cstddef>
#include <functional>


namespace  Geo
{

// Hash function
struct hashPairIntFunction
{
  size_t operator()(const std::pair<int , int> &x) const
  {
      int h1 = std::hash<int>{}(x.first);
      int h2 = std::hash<int>{}(x.second);
      return h1 ^ h2;
  }
};

inline float ToDeg(const float &a)
{
    return a*180.0f/float(M_PI);
}

}

#endif
