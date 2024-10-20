#ifndef T_JUNCTION
#define T_JUNCTION

namespace Geo {

//     VT
//      |
//      |
// V0---V---V1
struct TJunction 
{
  int V,V0,V1,VT;

  TJunction(){
    V=-1;
    V0=-1;
    V1=-1;
    VT=-1;
  }
  
  TJunction(int _V,int _V0,int _V1,int _VT)
  {
    V=_V;
    V0=_V0;
    V1=_V1;
    VT=_VT;
  }
};

} // namespace Geo
#endif
