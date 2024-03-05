#ifndef SMOOTH_MESH
#define SMOOTH_MESH

#include <cmath>
#include <vector>
#include "point3.h"
#include <solve/least_squares_lin_constraint_solver.h>
#include <connectivity.h>
#include <set>

namespace Geo
{


template <class ScalarType>
void SmoothLaplacian(std::vector<Geo::Point3<ScalarType> > &VertPos,
                     const std::vector<std::vector<int> > &Faces,int NumIte,
                     const std::vector<int> &FixedVert=std::vector<int>(),
                     const std::vector<std::pair<int,int> > &Features=std::vector<std::pair<int,int> >(),
                     bool line_border=true,ScalarType Damp=0.5)
{

    std::set<int> FixedVertSet(FixedVert.begin(),FixedVert.end());
    std::vector<bool> FeatureVert(VertPos.size(),false);
    for (size_t i=0;i<Features.size();i++)
    {
        FeatureVert[Features[i].first]=true;
        FeatureVert[Features[i].second]=true;
    }

    std::vector<bool> IsB;
    Geo::getBorderVert(Faces,IsB);
    for (int s=0;s<NumIte;s++)
    {
        std::vector<Geo::Point3<ScalarType> > NewVertPos(VertPos.size(),
                                                         Geo::Point3<ScalarType>(0,0,0));
        std::vector<int> SumNum(VertPos.size(),0);

        //cumulate features
        for (size_t i=0;i<Features.size();i++)
        {
            int IndexV0=Features[i].first;
            int IndexV1=Features[i].second;
            Geo::Point3<ScalarType> P0=VertPos[IndexV0];
            Geo::Point3<ScalarType> P1=VertPos[IndexV1];
            NewVertPos[IndexV1]+=P0;
            SumNum[IndexV1]++;
            NewVertPos[IndexV0]+=P1;
            SumNum[IndexV0]++;
        }

        for (size_t i=0;i<Faces.size();i++)
            for (size_t j=0;j<Faces[i].size();j++)
            {
                int IndexV0=Faces[i][j];
                int IndexV1=Faces[i][(j+1)%Faces[i].size()];
                Geo::Point3<ScalarType> P0=VertPos[IndexV0];
                Geo::Point3<ScalarType> P1=VertPos[IndexV1];

                if (line_border)
                {
                    if ((!IsB[IndexV0])&&(!FeatureVert[IndexV0]))
                    {
                        NewVertPos[IndexV0]+=P1;
                        SumNum[IndexV0]++;
                    }
                    if ((!IsB[IndexV1])&&(!FeatureVert[IndexV1]))
                    {
                        NewVertPos[IndexV1]+=P0;
                        SumNum[IndexV1]++;
                    }
                }
                else
                {
                    if (!FeatureVert[IndexV0])
                    {
                        NewVertPos[IndexV0]+=P1;
                        SumNum[IndexV0]++;
                    }
                    if (!FeatureVert[IndexV1])
                    {
                        NewVertPos[IndexV1]+=P0;
                        SumNum[IndexV1]++;
                    }
                }
            }
        if (line_border)
        {
            std::vector<std::pair<int,int> > BorderE;
            Geo::getBorderEdges(Faces,BorderE);
            for (size_t i=0;i<BorderE.size();i++)
            {
                int IndexV0=BorderE[i].first;
                int IndexV1=BorderE[i].second;
                Geo::Point3<ScalarType> P0=VertPos[IndexV0];
                Geo::Point3<ScalarType> P1=VertPos[IndexV1];

                if (!FeatureVert[IndexV0])
                {
                    NewVertPos[IndexV0]+=P1;
                    SumNum[IndexV0]++;
                }
                if (!FeatureVert[IndexV1])
                {
                    NewVertPos[IndexV1]+=P0;
                    SumNum[IndexV1]++;
                }
            }
        }
        for (size_t i=0;i<NewVertPos.size();i++)
        {
            if (FixedVertSet.count(i)>0)continue;
            if (SumNum[i]<2)continue;
            NewVertPos[i]/=(ScalarType)SumNum[i];
            VertPos[i]=NewVertPos[i]*(1-Damp)+VertPos[i]*Damp;
        }
    }
}

template <class ScalarType>
void SmoothLaplacianImplicit(std::vector<Geo::Point3<ScalarType> > &VertPos,
                             const std::vector<std::vector<int> > &Faces,
                             const std::vector<int> &FixedVert)
{
    ScalarType WeightFixed=1000000;
    ScalarType WeightLap=0.00001;

    std::vector<Geo::LinearConstraint<ScalarType> > LinConstr;
    std::set<int> FixedVertSet(FixedVert.begin(),FixedVert.end());
    //std::cout<<"there are "<<FixedVertSet.size()<<" fixed"<<std::endl;
    //add fixed velue
    for (size_t i=0;i<VertPos.size();i++)
    {
        if (FixedVertSet.count((int)i)>0)
        {
            for (size_t c=0;c<3;c++)
            {
                Geo::LinearConstraint<ScalarType> Constr;
                Constr.RHSide=VertPos[i].V(c);
                Constr.Weight=WeightFixed;
                Geo::VarTerm<ScalarType> Var(i,c,1);
                Constr.LHSide.push_back(Var);
                LinConstr.push_back(Constr);
            }
        }
    }

    //    //add laplacian along borders
    //    std::vector<std::pair<int,int> > BorderE;
    //    Geo::getBorderEdges(Faces,BorderE);

    //    for (size_t i=0;i<BorderE.size();i++)
    //    {
    //        int IndexV0=BorderE[i].first;
    //        int IndexV1=BorderE[i].second;

    //        for (size_t c=0;c<3;c++)
    //        {
    //            Geo::LinearConstraint<ScalarType> Constr;
    //            Constr.Weight=WeightBorder;
    //            Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV0,c,0.5));
    //            Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV1,c,-0.5));
    //            Constr.RHSide=0;
    //            LinConstr.push_back(Constr);
    //            //std::cout<<"ADDDED"<<std::endl;
    //        }
    //    }

    //add laplacian value
    std::vector<std::vector<int> > KernelV;
    Geo::FindVertexStar(Faces,VertPos.size(),KernelV);

    //    std::vector<bool> IsB;
    //    getBorderVert(Faces,IsB);

    for (size_t i=0;i<KernelV.size();i++)
    {
        int IndexV0=i;
        //        bool IsBV0=IsB[IndexV0];
        //        if (IsBV0)continue;
        for (size_t c=0;c<3;c++)
        {
            Geo::LinearConstraint<ScalarType> Constr;
            Constr.Weight=WeightLap;
            Constr.RHSide=0;
            Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV0,c,-1));
            ScalarType OtherW=(ScalarType)1/(ScalarType)KernelV[i].size();

            for (size_t j=0;j<KernelV[i].size();j++)
            {
                int IndexV1=KernelV[i][j];
                assert(IndexV0!=IndexV1);

                //                bool IsBV1=IsB[IndexV1];
                //                if (IsBV0 && (!IsBV1))continue;
                //                if (IsBV0 && IsBV1)
                //                    Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV1,c,0.5));
                //                else
                Constr.LHSide.push_back(Geo::VarTerm<ScalarType>(IndexV1,c,OtherW));
            }
            LinConstr.push_back(Constr);
        }
    }

    Geo::LinearConstraintLeastSQSolver<ScalarType> LinSqSolver;
    LinSqSolver.Init(LinConstr);
    std::vector<ScalarType> Solution;
    LinSqSolver.Solve(Solution);

    for (size_t i=0;i<VertPos.size();i++)
    {
        int IndexX=(i*3);
        int IndexY=IndexX+1;
        int IndexZ=IndexX+2;
        ScalarType ValX=Solution[IndexX];
        ScalarType ValY=Solution[IndexY];
        ScalarType ValZ=Solution[IndexZ];
        VertPos[i]=Geo::Point3<ScalarType >(ValX,ValY,ValZ);
    }
}


}
#endif
