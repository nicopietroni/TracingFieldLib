#ifndef POLYLINE_3_H
#define POLYLINE_3_H

#include <assert.h>
#include "point3.h"

namespace Geo
{

    template <class ScalarType>
    class Polyline3
    {
        std::vector<Geo::Point3<ScalarType>> PolyPos;

        void Add(const Geo::Point3<ScalarType> &P)
        {
            PolyPos.push_back(P);
        }
        
        // void SmoothStep()
        // {
        //     std::vector<Geo::Point3<ScalarType>> NewPolyPos;
        //     std::vector<ScalarType> SampleW;
        //     for (size_t i=1;i<PolyPos.size()-1;i++)
        //     {

        //     }
        // }

        void SmoothBilaplacian()
        {
            // for (size_t i=0;i<BorderVIdx.size();i++)
            //     for (size_t j=0;j<BorderVIdx[i].size();j++)
            //     {
            //         int SizeB=BorderVIdx[i].size();
            //         int Ind0=(j+(SizeB-2))%SizeB;
            //         int Ind1=(j+(SizeB-1))%SizeB;
            //         int Ind2=j;
            //         int Ind3=(j+1)%SizeB;
            //         int Ind4=(j+2)%SizeB;

            //         int VInd0=BorderVIdx[i][Ind0];
            //         int VInd1=BorderVIdx[i][Ind1];
            //         int VInd2=BorderVIdx[i][Ind2];
            //         int VInd3=BorderVIdx[i][Ind3];
            //         int VInd4=BorderVIdx[i][Ind4];

            //         //Geo::Point3<ScalarType> Pos0=TargetVertPos[VInd1];
            //         Geo::Point3<ScalarType> Pos1=TargetVertPos[VInd2];
            //         //Geo::Point3<ScalarType> Pos2=TargetVertPos[VInd3];
            //         //if (FixedPosSet.count(Pos0)>0)continue;
            //         if (FixedPosSet.count(Pos1)>0)continue;
            //         //if (FixedPosSet.count(Pos2)>0)continue;

            //         for (size_t c=0;c<3;c++)
            //         {
            //             Geo::LinearConstraint<ScalarType> Constr;
            //             Constr.Weight=WeightBorderSmooth;

            //             ScalarType TotW=0;
            //             Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(VInd0,c,0.1));
            //             Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(VInd1,c,0.4));
            //             Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(VInd2,c,-1));
            //             Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(VInd3,c,0.4));
            //             Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(VInd4,c,0.1));
            //             Constr.RHSide=0;
            //             LinConstr.push_back(Constr);
            //             numAdd++;
            //         }
            //     }
        }

        // ScalarType Lenght()
        // {
        //     ScalarType l = 0;
        //     for (size_t i = 0; i < PolyPos.size() - 1; i++)
        //     {
        //         Geo::Point2<ScalarType> P0 = PolyPos[i];
        //         Geo::Point2<ScalarType> P1 = PolyPos[i + 1];
        //         l += (P0 - P1).Norm();
        //     }
        //     return l;
        // }

        // void ParametrizeLenght(std::vector<ScalarType> &UVal)
        // {
        //     UVal.clear();
        //     ScalarType tot_l = Lenght();
        //     UVal.push_back(0);
        //     ScalarType curr_l = 0;
        //     for (size_t i = 0; i < PolyPos.size() - 1; i++)
        //     {
        //         Geo::Point2<ScalarType> P0 = PolyPos[i];
        //         Geo::Point2<ScalarType> P1 = PolyPos[i + 1];
        //         curr_l += (P0 - P1).Norm();
        //         UVal.push_back(curr_l / tot_l);
        //     }
        //     UVal.back()=1;
        // }

        // bool Interpolate(const ScalarType &AlphaVal,
        //                  int &IndexSeg,
        //                  ScalarType &SubAlpha)
        // {
        //     assert(AlphaVal >= 0);
        //     assert(AlphaVal <= 1);
        //     assert(Pos.size() >= 2);
        //     if (AlphaVal == 0)
        //     {
        //         IndexSeg = 0;
        //         SubAlpha = 0;
        //         return;
        //     }
        //     if (AlphaVal == 1)
        //     {
        //         IndexSeg = PolyPos.size() - 1;
        //         SubAlpha = 1;
        //         return;
        //     }
        //     std::vector<ScalarType> UVal;
        //     ParametrizeLenght(UVal);

        //     //then search the right one
        //     for (size_t i = 0; i < UVal.size() - 1; i++)
        //     {
        //         ScalarType U0=UVal[i];
        //         ScalarType U1=UVal[i+1];
        //         if ((AlphaVal>=U0)&&(AlphaVal<=U1))
        //         {
        //             IndexSeg=i;
        //             ScalarType Int0=AlphaVal-U0;
        //             ScalarType Int1=U1-AlphaVal;
        //             assert(Int0>=0);
        //             assert(Int1>=0);
        //             SubAlpha=1-(Int0/(Int0+Int1));
        //             return true;
        //         }
        //     }
        //     return false;
        // }
    };

}

#endif
