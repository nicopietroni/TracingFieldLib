#ifndef CONNECTIVITY
#define CONNECTIVITY

//#include "point3.h"
#include <Eigen/Sparse>
#include <set>
//#include "mesh_functions.h"

namespace Geo
{

//void ExtendStar( std::vector<std::vector<int> > &Kernel, int step)
//{
//    std::vector<std::vector<int>  > NewVertexKernels;
//    NewVertexKernels.resize(Kernel.size());
//    for (size_t i=0;i<Kernel.size();i++)
//    {
//        NewVertexKernels[i].push_back(i);
//        for (size_t j=0;j<Kernel[i].size();j++)
//        {
//            int NeighV=Kernel[i][j];
//            NewVertexKernels[i].push_back(NeighV);
//            NewVertexKernels[i].insert(NewVertexKernels[i].end(),
//                                       Kernel[NeighV].begin(),
//                                       Kernel[NeighV].end());
//        }
//    }
//    Kernel=NewVertexKernels;
//    for (size_t i=0;i<Kernel.size();i++)
//    {
//        std::sort(Kernel[i].begin(),Kernel[i].end());
//        auto last = std::unique(Kernel[i].begin(),Kernel[i].end());
//        Kernel[i].erase(last, Kernel[i].end());
//    }
//}

void ExtendStar( std::vector<std::vector<int> > &Kernel, int steps)
{
    std::vector<std::vector<int>  > Kernel0=Kernel;
    //NewVertexKernels.resize(Kernel.size());
    for (int s=0;s<steps;s++)
    {
        std::vector<std::vector<int>  > VertexKernel=Kernel;
        for (size_t i=0;i<VertexKernel.size();i++)
        {
            for (size_t j=0;j<VertexKernel[i].size();j++)
            {
                int NeighV=VertexKernel[i][j];
                for (size_t k=0;k<Kernel0[NeighV].size();k++)
                {
                    if (Kernel0[NeighV][k]!=(int)i)
                        Kernel[i].push_back(Kernel0[NeighV][k]);
                }
            }
        }
        //Kernel=NewVertexKernel;
        for (size_t i=0;i<Kernel.size();i++)
        {
            std::sort(Kernel[i].begin(),Kernel[i].end());
            auto last = std::unique(Kernel[i].begin(),Kernel[i].end());
            Kernel[i].erase(last, Kernel[i].end());
        }
    }
}

void FindVertexStar(const std::vector<std::vector<int> > &connections,
                    const size_t NumVert,
                    std::vector<std::vector<int> > &Kernel,
                    const int extend_step=0)
{
    Kernel.clear();
    Kernel.resize(NumVert);
    for (size_t i=0;i<connections.size();i++)
    {
        //int sizeF=connections[i].size();
        for (size_t j=0;j<connections[i].size();j++)
        {
            int VertIDX0=connections[i][j];
            for (size_t k=0;k<connections[i].size();k++)
            {
                if (j==k)continue;
                int VertIDX1=connections[i][k];
                assert(VertIDX0!=VertIDX1);
                Kernel[VertIDX0].push_back(VertIDX1);
                Kernel[VertIDX1].push_back(VertIDX0);
            }
        }
    }
    //make them unique
    for (size_t i=0;i<Kernel.size();i++)
    {
        std::sort(Kernel[i].begin(),Kernel[i].end());
        auto last = std::unique(Kernel[i].begin(),Kernel[i].end());
        Kernel[i].erase(last, Kernel[i].end());
    }

    //check symmetry
    for (size_t i=0;i<Kernel.size();i++)
    {
        int IndexV0=i;
        for (size_t j=0;j<Kernel[i].size();j++)
        {
            int IndexV1=Kernel[i][j];
            assert(std::find(Kernel[IndexV1].begin(), Kernel[IndexV1].end(), IndexV0) != Kernel[IndexV1].end());
        }
    }

    if (extend_step>0)
        ExtendStar(Kernel,extend_step);
}

int FindMaxVertIdx(const std::vector<std::vector<int> > &connections)
{
    int NumV=0;
    for (size_t i=0;i<connections.size();i++)
        for (size_t j=0;j<connections[i].size();j++)
            NumV=std::max(NumV,connections[i][j]);
    return NumV;
}

void FindVertexStar(const std::vector<std::vector<int> > &connections,
                    std::vector<std::vector<int> > &Kernel,
                    const int extend_step=0)
{
//    int NumV=0;
//    for (size_t i=0;i<connections.size();i++)
//        for (size_t j=0;j<connections[i].size();j++)
//            NumV=std::max(NumV,connections[i][j]);
    int NumV=FindMaxVertIdx(connections);
    NumV++;//num vertices is bigger index +1

    FindVertexStar(connections,NumV,Kernel,extend_step);
}

//void FindFacexStar(const std::vector<std::vector<int> > &Faces,
//                   const std::vector<std::vector<int> > &NextF,
//                   const std::vector<std::vector<int> > &NextE,
//                   std::vector<std::vector<int> > &Kernel,
//                   const int extend_step=0)
//{
//    //copy of nextf but remove borders
//    std::vector<std::vector<int> > NextFTemp;
//    NextFTemp.resize(Faces.size());
//    for (size_t i=0;i<NextF.size();i++)
//       for (size_t j=0;j<NextF[i].size();j++)
//            if (NextF[i][j]!=-1)NextFTemp[i].push_back(NextF[i][j]);

//    Kernel=NextFTemp;

//    for (size_t s=0;s<extend_step;s++)
//    {
//        for (size_t i=0;i<Kernel.size();i++)
//            for (size_t j=0;j<Kernel[i].size();j++)
//            {
//                int IndexF=Kernel[i][j];
//                assert(IndexF>=0);
//                assert(IndexF<NextFTemp.size());

//            }
//    }
////    FindVertexStar(connections,NumV,Kernel);

////    for (int i=0;i<extend_step;i++)
////        ExtendVertexStar(Kernel);
//}

void getFFAdjacency(const std::vector<std::vector<int> > &Faces,
                    std::vector<std::vector<int> > &NextF,
                    std::vector<std::vector<int> > &NextE)
{
    typedef typename std::pair<size_t,size_t> sizeTPair;

    NextF.clear();
    NextE.clear();
    NextF.resize(Faces.size());
    NextE.resize(Faces.size());

    std::map<sizeTPair, std::vector<sizeTPair> > EdgeVToFaceE;
    for (size_t i=0;i<Faces.size();i++)
    {
        size_t sizeF=Faces[i].size();
        NextF[i].resize(sizeF,-1);
        NextE[i].resize(sizeF,-1);
        for (size_t j=0;j<Faces[i].size();j++)
        {
            int indexV0=Faces[i][j];
            int indexV1=Faces[i][(j+1)%sizeF];
            sizeTPair Key((size_t)std::min(indexV0,indexV1),(size_t)std::max(indexV0,indexV1));
            EdgeVToFaceE[Key].push_back(sizeTPair(i,j));
        }
    }

    NextF.resize(Faces.size());
    NextE.resize(Faces.size());
    std::vector<std::pair<sizeTPair, std::vector<sizeTPair> > > EdgeVect(EdgeVToFaceE.begin(),EdgeVToFaceE.end());
    for (size_t i=0;i<EdgeVect.size();i++)
    {
        assert(EdgeVect[i].second.size()>=1);
        if (EdgeVect[i].second.size()==1)continue;//this is a border edge

        assert(EdgeVect[i].second.size()==2);

        int IndexF0=EdgeVect[i].second[0].first;
        int IndexE0=EdgeVect[i].second[0].second;

        int IndexF1=EdgeVect[i].second[1].first;
        int IndexE1=EdgeVect[i].second[1].second;

        assert(IndexF0!=IndexF1);

        assert(IndexE0<(int)NextF[IndexF0].size());
        assert(IndexE1<(int)NextF[IndexF1].size());

        NextF[IndexF0][IndexE0]=IndexF1;
        NextF[IndexF1][IndexE1]=IndexF0;

        NextE[IndexF0][IndexE0]=IndexE1;
        NextE[IndexF1][IndexE1]=IndexE0;
        //for (size_t j=0;j<EdgeVect.size();i++)

    }
}

void getVFAdjacency(const std::vector<std::vector<int> > &Faces,
                    std::vector<std::vector<int> > &VF,
                    std::vector<std::vector<int> > &VFIndex)
{
    int NumV=FindMaxVertIdx(Faces);

    VF.clear();
    VFIndex.clear();
    VF.resize(NumV+1);
    VFIndex.resize(NumV+1);

    for (size_t i=0;i<Faces.size();i++)
    {
        for (size_t j=0;j<Faces[i].size();j++)
        {
            int indexV=Faces[i][j];
            VF[indexV].push_back(i);
            VFIndex[indexV].push_back(j);
        }
    }
}

void getFFAdjacencyNonManifold(const std::vector<std::vector<int> > &Faces,
                               std::vector<std::vector<std::vector<int> > > &NextF,
                               std::vector<std::vector<std::vector<int> > > &NextE)
{
    typedef typename std::pair<size_t,size_t> sizeTPair;

    NextF.clear();
    NextE.clear();
    NextF.resize(Faces.size());
    NextE.resize(Faces.size());

    std::map<sizeTPair, std::vector<sizeTPair> > EdgeVToFaceE;
    for (size_t i=0;i<Faces.size();i++)
    {
        size_t sizeF=Faces[i].size();
        NextF[i].resize(sizeF);
        NextE[i].resize(sizeF);
        for (size_t j=0;j<Faces[i].size();j++)
        {
            int indexV0=Faces[i][j];
            int indexV1=Faces[i][(j+1)%sizeF];
            sizeTPair Key((size_t)std::min(indexV0,indexV1),(size_t)std::max(indexV0,indexV1));
            EdgeVToFaceE[Key].push_back(sizeTPair(i,j));
        }
    }

    std::vector<std::pair<sizeTPair, std::vector<sizeTPair> > > EdgeVect(EdgeVToFaceE.begin(),EdgeVToFaceE.end());
    for (size_t i=0;i<EdgeVect.size();i++)
    {
        assert(EdgeVect[i].second.size()>=1);
        if (EdgeVect[i].second.size()==1)//this is a border edge
        {
            int IndexF0=EdgeVect[i].second[0].first;
            int IndexE0=EdgeVect[i].second[0].second;
            NextF[IndexF0][IndexE0].push_back(-1);
            NextE[IndexF0][IndexE0].push_back(-1);
        }
        else
        {
            for(size_t j=0;j<EdgeVect[i].second.size();j++)
            {
                int IndexF0=EdgeVect[i].second[j].first;
                int IndexE0=EdgeVect[i].second[j].second;
                for(size_t k=0;k<EdgeVect[i].second.size();k++)
                {
                    if (j==k)continue;
                    int IndexF1=EdgeVect[i].second[k].first;
                    int IndexE1=EdgeVect[i].second[k].second;
                    assert(IndexF0!=IndexF1);
                    NextF[IndexF0][IndexE0].push_back(IndexF1);
                    NextE[IndexF0][IndexE0].push_back(IndexE1);
                }
            }
        }
    }
}

void FindFaceStar(const std::vector<std::vector<int> > &connections,
                  std::vector<std::vector<int> > &Kernel,
                  const int extend_step=0)
{
    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(connections, NextF,NextE);

    Kernel.resize(NextF.size());

    for (size_t i=0;i<NextF.size();i++)
        for (size_t j=0;j<NextF[i].size();j++)
            if (NextF[i][j]!=-1)Kernel[i].push_back(NextF[i][j]);

//    for (int i=0;i<extend_step;i++)
    ExtendStar(Kernel,extend_step);
}

void getBorderVert(const std::vector<std::vector<int> > &connections,
                   const std::vector<std::vector<int> > &NextF,
                   std::vector<size_t> &BorderV)
{
    BorderV.clear();
    for (size_t i=0;i<NextF.size();i++)
        for (size_t j=0;j<NextF[i].size();j++)
        {
            if (NextF[i][j]!=-1)continue;
            int numV=connections[i].size();
            int IndexV0=connections[i][j];
            int IndexV1=connections[i][(j+1)%numV];
            BorderV.push_back(IndexV0);
            BorderV.push_back(IndexV1);
        }
    std::sort(BorderV.begin(),BorderV.end());
    std::vector<size_t>::iterator last = std::unique(BorderV.begin(), BorderV.end());
    BorderV.erase(last, BorderV.end());
}

void getBorderVert(const std::vector<std::vector<int> > &connections,
                   std::vector<size_t> &BorderV)
{
    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(connections, NextF,NextE);
    getBorderVert(connections,NextF,BorderV);
}

void getBorderEdges(const std::vector<std::vector<int> > &Faces,
                    const std::vector<std::vector<int> > &NextF,
                    std::vector<std::pair<int,int> > &BorderE)
{
    BorderE.clear();
    for (size_t i=0;i<NextF.size();i++)
        for (size_t j=0;j<NextF[i].size();j++)
        {
            if (NextF[i][j]!=-1)continue;
            int numV=Faces[i].size();
            int IndexV0=Faces[i][j];
            int IndexV1=Faces[i][(j+1)%numV];
            std::pair<int,int> EdgeKey(std::min(IndexV0,IndexV1),
                                       std::max(IndexV0,IndexV1));
            BorderE.push_back(EdgeKey);
        }
    std::sort(BorderE.begin(),BorderE.end());
    std::vector<std::pair<int,int> >::iterator last = std::unique(BorderE.begin(), BorderE.end());
    BorderE.erase(last, BorderE.end());
}


void getBorderEdges(const std::vector<std::vector<int> > &Faces,
                    std::vector<std::pair<int,int> > &BorderE)
{
    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(Faces,NextF,NextE);
    getBorderEdges(Faces,NextF,BorderE);
}

void getBorderVert(const std::vector<std::vector<int> > &connections,
                   const std::vector<std::vector<int> > &NextF,
                   const int NumVert,
                   std::vector<bool> &IsB)
{
    std::vector<size_t> BorderV;
    getBorderVert(connections,NextF,BorderV);

    IsB.clear();
    IsB.resize(NumVert,false);
    for (size_t i=0;i<BorderV.size();i++)
    {
        size_t IndexV=BorderV[i];
        assert(IndexV<IsB.size());
        IsB[IndexV]=true;
    }
}


void getBorderVert(const std::vector<std::vector<int> > &connections,
                   std::vector<bool> &IsB)
{
    std::vector<std::vector<int> > NextF,NextE;
    Geo::getFFAdjacency(connections, NextF,NextE);
    int NumV=Geo::FindMaxVertIdx(connections);
    getBorderVert(connections,NextF,NumV+1,IsB);
}


}

#endif
