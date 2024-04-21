#ifndef POLYLINE_2_H
#define POLYLINE_2_H

#include <assert.h>
#include "point2.h"

namespace Geo
{

    // template <class ScalarType>
    // class Polyline2
    // {
    //     std::vector<Geo::Point2<ScalarType>> PolyPos;

    //     void Add(const Geo::Point2<ScalarType> &P)
    //     {
    //         PolyPos.push_back(P);
    //     }

    //     ScalarType Lenght()
    //     {
    //         ScalarType l = 0;
    //         for (size_t i = 0; i < PolyPos.size() - 1; i++)
    //         {
    //             Geo::Point2<ScalarType> P0 = PolyPos[i];
    //             Geo::Point2<ScalarType> P1 = PolyPos[i + 1];
    //             l += (P0 - P1).Norm();
    //         }
    //         return l;
    //     }

    //     void ParametrizeLenght(std::vector<ScalarType> &UVal)
    //     {
    //         UVal.clear();
    //         ScalarType tot_l = Lenght();
    //         UVal.push_back(0);
    //         ScalarType curr_l = 0;
    //         for (size_t i = 0; i < PolyPos.size() - 1; i++)
    //         {
    //             Geo::Point2<ScalarType> P0 = PolyPos[i];
    //             Geo::Point2<ScalarType> P1 = PolyPos[i + 1];
    //             curr_l += (P0 - P1).Norm();
    //             UVal.push_back(curr_l / tot_l);
    //         }
    //         UVal.back()=1;
    //     }

    //     bool Interpolate(const ScalarType &AlphaVal,
    //                      int &IndexSeg,
    //                      ScalarType &SubAlpha)
    //     {
    //         assert(AlphaVal >= 0);
    //         assert(AlphaVal <= 1);
    //         assert(Pos.size() >= 2);
    //         if (AlphaVal == 0)
    //         {
    //             IndexSeg = 0;
    //             SubAlpha = 0;
    //             return;
    //         }
    //         if (AlphaVal == 1)
    //         {
    //             IndexSeg = PolyPos.size() - 1;
    //             SubAlpha = 1;
    //             return;
    //         }
    //         std::vector<ScalarType> UVal;
    //         ParametrizeLenght(UVal);

    //         //then search the right one
    //         for (size_t i = 0; i < UVal.size() - 1; i++)
    //         {
    //             ScalarType U0=UVal[i];
    //             ScalarType U1=UVal[i+1];
    //             if ((AlphaVal>=U0)&&(AlphaVal<=U1))
    //             {
    //                 IndexSeg=i;
    //                 ScalarType Int0=AlphaVal-U0;
    //                 ScalarType Int1=U1-AlphaVal;
    //                 assert(Int0>=0);
    //                 assert(Int1>=0);
    //                 SubAlpha=1-(Int0/(Int0+Int1));
    //                 return true;
    //             }
    //         }
    //         return false;
    //     }
    // };

}

#endif
