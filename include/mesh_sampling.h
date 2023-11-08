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
    static void GetRandomSeq(size_t MaxIndex,std::vector<size_t> &Sequence,int RandSeed=0)
    {
        Sequence.clear();
        for (size_t i = 0; i <= MaxIndex; i++) Sequence.push_back(i);

        std::mt19937 gen(RandSeed);
        std::shuffle(Sequence.begin(), Sequence.end(), gen);
    }

public:

    static void SampleVert(const std::vector<std::vector<int> > &Connectivity,
                           const int NumVert,const int &RadiusStar,
                           std::vector<int> &Selected,int RandSeed=0)
    {
        Selected.clear();
        if (RadiusStar==0)
        {
            for (int i = 0; i < NumVert; i++) Selected.push_back(i);
            return;
        }
        //initialize the kernel
        std::vector<std::vector<int> > KernelVert;
        Geo::FindVertexStar(Connectivity,KernelVert,RadiusStar-1);

        std::vector<size_t> Sequence;
        std::vector<bool> Valid(NumVert,true);
        //randomly perturb the order
        GetRandomSeq(NumVert-1,Sequence,RandSeed);

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

    static void SampleFaces(const std::vector<std::vector<int> > &Connectivity,
                            const int &RadiusStar,std::vector<int> &Selected,
                            int RandSeed=0)
    {
        Selected.clear();
        if (RadiusStar==0)
        {
            for (int i = 0; i < Connectivity.size(); i++) Selected.push_back(i);
            return;
        }
        //initialize the kernel
        std::vector<std::vector<int> > KernelVert;
        Geo::FindFaceStar(Connectivity,KernelVert,RadiusStar-1);

        std::vector<size_t> Sequence;
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
