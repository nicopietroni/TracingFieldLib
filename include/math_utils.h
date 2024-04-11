#ifndef MATH_UTILS
#define MATH_UTILS

#include <cmath>
#include <utility>
#include <cstddef>
#include <functional>
//#include <math.h>
#define M_PI        3.14159265358979323846264338327950288   /* pi             */

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

template <class ScalarType>
inline ScalarType ToDeg(const ScalarType &a)
{
    return a*(ScalarType)180.0f/ScalarType(M_PI);
}


template <class ScalarType>
inline ScalarType normalizeAngle(ScalarType angle) 
{
    angle =(ScalarType) fmod((float)angle, (float) 2 * M_PI);
    if (angle < 0) {
        angle += 2 * M_PI;
    }
    return angle;
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

template <class ScalarType>
inline void Trim(std::vector<ScalarType> &Values,
                 const ScalarType &MaxVal)
{
    for (size_t i=0;i<Values.size();i++)
        Values[i]=std::min(Values[i],MaxVal);
}


}

#endif
