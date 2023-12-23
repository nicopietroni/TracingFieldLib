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

template <class ScalarType>
inline ScalarType MaxPercentile(const std::vector<ScalarType> &Values,
                                const ScalarType &Percentile)
{
    assert((Percentile>=0)&&(Percentile<=1));

    //then compute the max (considering the percentile)
    std::vector<ScalarType> ValuesSorted=Values;
    std::sort(ValuesSorted.begin(),ValuesSorted.end());
    int IndexMax=floor(((ScalarType)ValuesSorted.size()*((ScalarType)1-Percentile))+0.5);
    IndexMax=std::min(IndexMax,(int)ValuesSorted.size()-1);
    IndexMax=std::max(IndexMax,0);
    ScalarType TestVal=ValuesSorted[IndexMax];
    return TestVal;
}


}

#endif
