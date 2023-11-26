#ifndef MESH_SAMPLING
#define SAMPLE_FOR_VIS_FRAME_FIELD

//#include "point3.h"
//#include "polyhedra_mesh.h"
#include <algorithm>
#include <random>
#include <vector>
#include "connectivity.h"

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
                             const int &RadiusStar,std::vector<int> &To_prune)
    {
        //initialize the kernel
        std::vector<std::vector<int> > KernelVert;
        Geo::FindVertexStar(Connectivity,KernelVert,RadiusStar-1);

        std::vector<int> Selected;
        std::vector<bool> Valid(NumVert,true);
        for (size_t i=0;i<To_prune.size();i++)
        {
            size_t CurrV=To_prune[i];
            if (!Valid[CurrV])continue;

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

//    static void SampleBorder(const std::vector<std::vector<int> > &Connectivity,
//                           const int NumVert,const int &RadiusStar,
//                           std::vector<size_t> &Selected,int RandSeed=0)
//    {
//        //get all border
//        std::vector<std::vector<int> > NextF,NextE;
//        Geo::getFFAdjacency(Connectivity,NextF,NextE);
//        std::vector<size_t> BorderV;
//        Geo::getBorderVert(Connectivity,NextF,BorderV);
//        Selected=BorderV;

//        //shuffle
//        std::mt19937 gen(RandSeed);
//        std::shuffle(Selected.begin(), Selected.end(), gen);

//        //then prune
//        PruneSamples(NumVert,Connectivity,RadiusStar,Selected);
//    }

//    static void SampleVert(const std::vector<std::vector<int> > &Connectivity,
//                           const int NumVert,const int &RadiusStar,
//                           std::vector<int> &Selected,
//                           int RandSeed=0)
//    {
//        Selected.clear();
//        if (RadiusStar==0)
//        {
//            for (int i = 0; i < NumVert; i++) Selected.push_back(i);
//            return;
//        }
//        //initialize the kernel
//        std::vector<std::vector<int> > KernelVert;
//        Geo::FindVertexStar(Connectivity,KernelVert,RadiusStar-1);

//        std::vector<size_t> Sequence;
//        std::vector<bool> Valid(NumVert,true);

//        //randomly perturb the order
//        GetRandomSeq(NumVert-1,Sequence,RandSeed);

//        for (size_t i=0;i<Sequence.size();i++)
//        {
//            int CurrV=Sequence[i];
//            if (!Valid[CurrV])continue;

//            Selected.push_back(CurrV);
//            for (size_t j=0;j<KernelVert[CurrV].size();j++)
//            {
//                int NeighV=KernelVert[CurrV][j];
//                Valid[NeighV]=false;
//            }
//        }
//    }

    static void SampleVert(const std::vector<std::vector<int> > &Connectivity,
                           const int NumVert,const int &RadiusStar,
                           std::vector<int> &Selected,bool priority_border,
                           int RandSeed=0)
    {
        Selected.clear();
        if (RadiusStar==0)
        {
            for (int i = 0; i < NumVert; i++) Selected.push_back(i);
            return;
        }

        if (priority_border)
        {
            std::vector<std::vector<int> > NextF,NextE;
            Geo::getFFAdjacency(Connectivity,NextF,NextE);
            std::vector<bool> IsB;
            Geo::getBorderVert(Connectivity,NextF,NumVert,IsB);
            assert(IsB.size()==NumVert);

            std::vector<int> NonBordV,BorderV;

            for (size_t i=0;i<IsB.size();i++)
            {
                if (IsB[i])
                    BorderV.push_back((int)i);
                else
                    NonBordV.push_back((int)i);
            }
            std::mt19937 gen(RandSeed);
            std::shuffle(BorderV.begin(), BorderV.end(), gen);
            std::shuffle(NonBordV.begin(), NonBordV.end(), gen);
            //first border then others
            Selected=BorderV;
            Selected.insert(Selected.end(),NonBordV.begin(),NonBordV.end());
        }
        else
        {
            GetRandomSeq(NumVert-1,Selected,RandSeed);
        }
        PruneSamples(NumVert,Connectivity,RadiusStar,Selected);
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

};

}

#endif
