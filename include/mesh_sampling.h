#ifndef MESH_SAMPLING
#define MESH_SAMPLING

//#include "point3.h"
//#include "polyhedra_mesh.h"
#include <algorithm>
#include <random>
#include <vector>
#include "connectivity.h"
#include "point2.h"
#include <pos3.h>
#include <Space_Query/point_grid_2D.h>

namespace Geo{

template <class ScalarType>
class MeshSampling
{
    static void GetRandomSeq(size_t MaxIndex,std::vector<int> &Sequence,int RandSeed=0)
    {
        Sequence.clear();
        for (size_t i = 0; i <= MaxIndex; i++) Sequence.push_back(i);

        std::mt19937 gen(RandSeed);
        std::shuffle(Sequence.begin(), Sequence.end(), gen);
    }

    static void PruneSamples(const int NumVert,
                             const std::vector<std::vector<int> > &Connectivity,
                             const int &RadiusStar,
                             const std::vector<int> &Constraints,
                             std::vector<int> &To_prune)
    {
        //initialize the kernel
        std::vector<std::vector<int> > KernelVert;
        Geo::FindVertexStar(Connectivity,KernelVert,RadiusStar-1);

        std::set<int> ConstraintsSet(Constraints.begin(),Constraints.end());

        std::vector<int> Selected;
        std::vector<bool> Valid(NumVert,true);

        //first prune considering the border
        for (size_t i=0;i<Constraints.size();i++)
        {
            //add those by default
            size_t CurrV=Constraints[i];
            Selected.push_back(CurrV);

            //invalidate the rest
            for (size_t j=0;j<KernelVert[CurrV].size();j++)
            {
                int NeighV=KernelVert[CurrV][j];
                Valid[NeighV]=false;
            }
        }

        for (size_t i=0;i<To_prune.size();i++)
        {
            size_t CurrV=To_prune[i];

            //this already inserted and pruned
            if (ConstraintsSet.count(CurrV)>0)
                continue;

            //other one non valid, go on
            if (!Valid[CurrV])continue;

            //add those as selected
            Selected.push_back(CurrV);
            for (size_t j=0;j<KernelVert[CurrV].size();j++)
            {
                int NeighV=KernelVert[CurrV][j];
                Valid[NeighV]=false;
            }
        }
        To_prune=Selected;
    }


public:

    static void SampleVert(const std::vector<std::vector<int> > &Connectivity,
                           const int NumVert,const int &RadiusStar,
                           std::vector<int> &Selected,
                           const std::vector<int> &Constraints=std::vector<int>(),
                           int RandSeed=0,std::vector<int> *SampleFrom=NULL)
    {
        Selected.clear();
        if (RadiusStar==0)
        {
            if (SampleFrom!=NULL)
            {
                for (int i = 0; i < (*SampleFrom).size(); i++)
                    Selected.push_back((*SampleFrom)[i]);
            }
            else
                for (int i = 0; i < NumVert; i++)
                    Selected.push_back(i);

            return;
        }

        if (SampleFrom!=NULL)
        {
            std::mt19937 gen(RandSeed);
            std::shuffle((*SampleFrom).begin(),(*SampleFrom).end(), gen);
            Selected=(*SampleFrom);
            PruneSamples(NumVert,Connectivity,RadiusStar,Constraints,Selected);
        }
        else
        {
            GetRandomSeq(NumVert-1,Selected,RandSeed);
            PruneSamples(NumVert,Connectivity,RadiusStar,Constraints,Selected);
       }
    }

    static void SampleFaces(const std::vector<std::vector<int> > &Connectivity,
                            const int &RadiusStar,std::vector<int> &Selected,
                            int RandSeed=0)
    {
        Selected.clear();
        if (RadiusStar==0)
        {
            for (int i = 0; i < (int)Connectivity.size(); i++) Selected.push_back(i);
            return;
        }
        //initialize the kernel

        std::vector<std::vector<int> > KernelVert;
        Geo::FindFaceStar(Connectivity,KernelVert,RadiusStar-1);

        std::vector<int> Sequence;
        std::vector<bool> Valid(Connectivity.size(),true);
        //randomly perturb the order
        GetRandomSeq(Connectivity.size()-1,Sequence,RandSeed);

        for (size_t i=0;i<Sequence.size();i++)
        {
            int CurrV=Sequence[i];
            if (!Valid[CurrV])continue;

            Selected.push_back(CurrV);
            for (size_t j=0;j<KernelVert[CurrV].size();j++)
            {
                int NeighV=KernelVert[CurrV][j];
                Valid[NeighV]=false;
            }
        }
    }

    static void SampleRadiusUV(const std::vector<Geo::Point2<ScalarType> > &Samples,
                               const ScalarType &RadiusStar,std::vector<int> &Selected,
                                int RandSeed=0)
    {
        Geo::PointSpatialIndex2D<ScalarType> PointIdx;
        PointIdx.Init(Samples);

        Selected.clear();

        std::vector<bool> Valid(Samples.size(),true);

        //randomly perturb the order
        std::vector<int> Sequence;
        GetRandomSeq(Samples.size()-1,Sequence,RandSeed);

        for (size_t i=0;i<Sequence.size();i++)
        {
            int CurrV=Sequence[i];
            if (!Valid[CurrV])continue;

            Selected.push_back(CurrV);
            std::vector<int> PosIn;
            PointIdx.PointsInRadius(Samples,Samples[CurrV],RadiusStar,PosIn);

            for (size_t j=0;j<PosIn.size();j++)
            {
                int CurrV1=PosIn[j];
                if (CurrV1!=CurrV)
                    Valid[CurrV1]=false;
            }
        }
    }

    static void SampleBordersRegularly(const std::vector<std::vector<int> > &Faces,
                                       const std::vector<std::vector<int> > &NextF,
                                       const std::vector<std::vector<int> > &NextE,
                                       const int StepSize,const int MinSamplePerBorder,
                                       std::vector<int> &Selected,const std::vector<int> *SampleFrom=NULL)

    {
        std::set<int> SampleFromSet;
        if (SampleFrom!=NULL)
            SampleFromSet=std::set<int>((*SampleFrom).begin(),SampleFrom->end());

        std::vector<std::vector<Pos3 > > BorderPos;
        getBorderSequnces(Faces,NextF,NextE,BorderPos);
        for (size_t i=0;i<BorderPos.size();i++)
        {
            int CurrStepSize=StepSize;
            int Subd=floor((ScalarType)BorderPos[i].size()/(ScalarType)StepSize);
            if (Subd<MinSamplePerBorder)
            {
                CurrStepSize=floor((ScalarType)BorderPos[i].size()/(ScalarType)MinSamplePerBorder);
            }
            CurrStepSize=std::max((int)1,CurrStepSize);

            for (size_t j=0;j<BorderPos[i].size();j+=CurrStepSize)
            {
                if (SampleFrom!=NULL)
                {
                    if (SampleFromSet.count(BorderPos[i][j].V)>0)
                        Selected.push_back(BorderPos[i][j].V);
                }
                else
                    Selected.push_back(BorderPos[i][j].V);
            }
        }
    }

    static void SampleWithBorderRegular(const int NumVert,
                                        const std::vector<std::vector<int> > &Faces,
                                        const std::vector<std::vector<int> > &NextF,
                                        const std::vector<std::vector<int> > &NextE,
                                        const int StepSize, const int MinSamplePerBorder,
                                        std::vector<int> &Selected,
                                        std::vector<int> *SampleFrom=NULL,
                                        int RandSeed=0)
    {
        std::vector<int> BorderConstraints;
        SampleBordersRegularly(Faces,NextF,NextE,StepSize,MinSamplePerBorder,BorderConstraints,SampleFrom);
        SampleVert(Faces,NumVert,StepSize,Selected,BorderConstraints,RandSeed,SampleFrom);
    }

    static void SampleWithBorderRegular(const int NumVert,
                                        const std::vector<std::vector<int> > &Faces,
                                        const int StepSize, const int MinSamplePerBorder,
                                        std::vector<int> &Selected,
                                        std::vector<int> *SampleFrom=NULL,
                                        int RandSeed=0)
    {
        std::vector<std::vector<int> > NextF,NextE;
        Geo::getFFAdjacency(Faces,NextF,NextE);
        SampleWithBorderRegular(NumVert,Faces,NextF,NextE,StepSize,MinSamplePerBorder,Selected,SampleFrom,RandSeed);
    }

    static int EsteemSamplingRadius(const int NumVer,const int TargetSamples)
    {
        ScalarType AreaC=(ScalarType)NumVer/(ScalarType)TargetSamples;
        int RadiusStar=floor(sqrt(AreaC)+0.5);
        return RadiusStar;
    }
};

}

#endif
