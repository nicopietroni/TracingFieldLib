#ifndef COLOR_FUNCTIONS
#define COLOR_FUNCTIONS

#include <vector>
#include "point3.h"

namespace Geo
{

//input: ratio is between 0.0 to 1.0
//output: rgb color
template <class ScalarType>
Point3<ScalarType> InterpColor(const Point3<ScalarType> &Col0,
                               const Point3<ScalarType> &Col1,
                               const ScalarType &ratio)
{
    assert(ratio>=0);
    assert(ratio<=1);

    Point3<ScalarType> RetCol=Col1*(ratio)+Col0*(1-ratio);
    return RetCol;
}

template <class ScalarType>
Point3<ScalarType> HSVColor( float h, float s, float v)
{
    float r,g,b;
    Point3<ScalarType> RetCol;
    if(s==0.0){	// gray color
        r = g = b = v;
        RetCol.X=(r);
        RetCol.Y=(g);
        RetCol.Z=(b);
        return RetCol;
    }
    float dummy;
    h = modff(h,&dummy);
    if(h==1.0) h = 0.0;

    int i   = int( floor(h*6.0f) );
    float f = float(h*6.0f- floor(h*6.0f));

    float p = v*(1.0f-s);
    float q = v*(1.0f-s*f);
    float t = v*(1.0f-s*(1.0f-f));

    switch(i)
    {
    case 0: r=v; g=t; b=p; break;
    case 1: r=q; g=v; b=p; break;
    case 2: r=p; g=v; b=t; break;
    case 3: r=p; g=q; b=v; break;
    case 4: r=t; g=p; b=v; break;
    case 5: r=v; g=p; b=q; break;
    default: r=0;g=0;b=0; assert(0);break;
    }
    RetCol.X=(r);
    RetCol.Y=(g);
    RetCol.Z=(b);
    return RetCol;
}

template <class ScalarType>
Point3<ScalarType> ScatterColor(int range, int value,float Sat=.3f,float Val=.9f)
{
    int b, k, m=range;
    int r =range;

    for (b=0, k=1; k<range; k<<=1)
        if (value<<1>=m) {
            if (b==0) r = k;
            b += k;
            value -= (m+1)>>1;
            m >>= 1;
        }
        else m = (m+1)>>1;

    if (r>range-b) r = range-b;

    //TRACE("Scatter range 0..%i, in %i out %i\n",n,a,b);
    //SetHSVColor(float(b)/float(range),Sat,Val);
    return (HSVColor<ScalarType>(float(b)/float(range),Sat,Val));
}

template <class ScalarType>
inline Geo::Point3<ScalarType> Blue()
{
    return Geo::Point3<ScalarType>(0,0,1);
}

template <class ScalarType>
inline Geo::Point3<ScalarType> Green()
{
    return Geo::Point3<ScalarType>(0,1,0);
}

template <class ScalarType>
inline Geo::Point3<ScalarType> Cyan()
{
    return Geo::Point3<ScalarType>(0,1,1);
}

template <class ScalarType>
inline Geo::Point3<ScalarType> Yellow()
{
    return Geo::Point3<ScalarType>(1,1,0);
}

template <class ScalarType>
inline Geo::Point3<ScalarType> Red()
{
    return Geo::Point3<ScalarType>(1,0,0);
}

template <class ScalarType>
inline Geo::Point3<ScalarType> ColorRGBRamp(const ScalarType minf,const ScalarType  maxf ,ScalarType v )
{
    if(minf>maxf) { return(ColorRGBRamp(maxf,minf,maxf+(minf-v))); }

    ScalarType step=(maxf-minf)/4;
    if(v <  minf ) { return Red<ScalarType>(); }
    v-=minf;
    if(v<step)
    {

        return (InterpColor(Red<ScalarType>(),Yellow<ScalarType>(),v/step));
    }
    v-=step;
    if(v<step)
    {
        return (InterpColor(Yellow<ScalarType>(),Green<ScalarType>(),v/step));
    }
    v-=step;
    if(v<step)
    {
        return (InterpColor(Green<ScalarType>(),Cyan<ScalarType>(),v/step));
    }
    v-=step;
    if(v<step)
    {
        return (InterpColor(Cyan<ScalarType>(),Blue<ScalarType>(),v/step));
    }
    return Blue<ScalarType>();
}

template <class ScalarType>
void GetColorByScalar(const std::vector<ScalarType> &ScalarVal,
                      ScalarType MinVal,ScalarType MaxVal,
                      std::vector<Geo::Point3<ScalarType> > &Color)
{
    assert(MinVal!=MaxVal);
    Color.resize(ScalarVal.size(),Geo::Point3<ScalarType> (1,1,1));
    for (size_t i=0;i<Color.size();i++)
    {
        ScalarType Val=ScalarVal[i];
        Color[i]=ColorRGBRamp(MinVal,MaxVal ,Val );
    }
}

template <class ScalarType>
void GetColorByScalar(const std::vector<ScalarType> &ScalarVal,
                      std::vector<Geo::Point3<ScalarType> > &Color)
{
    // Find the minimum element
    typename std::vector<ScalarType>::const_iterator minElement = std::min_element(ScalarVal.begin(), ScalarVal.end());
    // Find the maximum element
    typename std::vector<ScalarType>::const_iterator maxElement = std::max_element(ScalarVal.begin(), ScalarVal.end());
    ScalarType MinVal=(*minElement);
    ScalarType MaxVal=(*maxElement);
    GetColorByScalar(ScalarVal,MinVal,MaxVal,Color);

}
}

#endif
