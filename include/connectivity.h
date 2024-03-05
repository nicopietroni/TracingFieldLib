#ifndef CONNECTIVITY
#define CONNECTIVITY

//#include "point3.h"
#include <Eigen/Sparse>
#include <iostream>
#include <set>
//#include "mesh_functions.h"

namespace Geo
{


int WhichVert(const std::vector<int> &Face,
              const int &VertIndex)
{
    for (size_t i=0;i<Face.size();i++)
        if (Face[i]==VertIndex)return i;
    return -1;
}

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


bool CheckFFAdjacencyManifold(const std::vector<std::vector<int> > &Faces,
                              std::vector<std::vector<int> > &NextF,
                              std::vector<std::vector<int> > &NextE)
{
    //check entries to faces
    for (size_t i=0;i<NextF.size();i++)
        for (size_t j=0;j<NextF[i].size();j++)
        {
            int IndexF0=i;
            int IndexE0=j;
            int IndexF1=NextF[i][j];
            int IndexE1=NextE[i][j];
            if (IndexF1==-1)
            {
                if (IndexE1!=-1)return false;
                continue;
            }

            //check vertices
            int SizeF0=Faces[IndexF0].size();
            int IndexV0_F0=Faces[IndexF0][IndexE0];
            int IndexV1_F0=Faces[IndexF0][(IndexE0+1)%SizeF0];

            int SizeF1=Faces[IndexF1].size();
            int IndexV0_F1=Faces[IndexF1][IndexE1];
            int IndexV1_F1=Faces[IndexF1][(IndexE1+1)%SizeF1];

            if(IndexV0_F0!=IndexV1_F1)return false;
            if(IndexV1_F0!=IndexV0_F1)return false;
        }

    //faces to entry
    std::vector<std::vector<int> > NextFTest;
    std::vector<std::vector<int> > NextETest;
    getFFAdjacency(Faces,NextFTest,NextETest);
    if(NextF!=NextFTest)return false;
    if(NextE!=NextETest)return false;

    return true;
}

bool HasVFAdj(const std::vector<std::vector<int> > &VF,
              const int &VIndex,const int &FIndex)
{
    assert(VIndex>=0);
    assert(VIndex<(int)VF.size());
    std::vector<int>::const_iterator it = std::find(VF[VIndex].begin(),
                                                    VF[VIndex].end(),
                                                    FIndex);

    return(it != VF[VIndex].end());
}

void VFDetach(std::vector<std::vector<int> > &VF,
              std::vector<std::vector<int> > &VFIndex,
              const int &VIndex,
              const int &FIndex)
{
    assert(VFIndex.size()==VF.size());
    assert(VIndex<(int)VF.size());
    assert(VIndex<(int)VFIndex.size());
    std::vector<int> newVF;
    std::vector<int> newVFIdx;
    for (size_t i=0;i<VF[VIndex].size();i++)
    {
        if (VF[VIndex][i]!=FIndex)
        {
            newVF.push_back(VF[VIndex][i]);
            newVFIdx.push_back(VFIndex[VIndex][i]);
        }
    }
    //    std::cout<<"Size Before:"<<VF[VIndex].size()<<std::endl;
    //    std::cout<<"Size After:"<<newVF.size()<<std::endl;
    assert(newVF.size()==VF[VIndex].size()-1);
    assert(newVFIdx.size()==VFIndex[VIndex].size()-1);
    VF[VIndex]=newVF;
    VFIndex[VIndex]=newVFIdx;
}

void VFAttach(const std::vector<std::vector<int> > &Faces,
              std::vector<std::vector<int> > &VF,
              std::vector<std::vector<int> > &VFIndex,
              const int &VIndex,
              const int &FIndex)
{
    assert(VFIndex.size()==VF.size());
    assert(VIndex<(int)VF.size());
    assert(VIndex<(int)VFIndex.size());
    assert(std::find(VF[VIndex].begin(),VF[VIndex].end(), FIndex)==VF[VIndex].end());
    VF[VIndex].push_back(FIndex);
    int LocalIdx=WhichVert(Faces[FIndex],VIndex);
    assert(LocalIdx>=0);
    assert(LocalIdx<3);
    VFIndex[VIndex].push_back(LocalIdx);
}

void UpdateVF(const std::vector<std::vector<int> > &Faces,
              const std::vector<int> &To_UpdF,
              std::vector<std::vector<int> > &VF,
              std::vector<std::vector<int> > &VFIndex)
{
    //get all vert involved
    std::set<int> To_UpdVSet;
    for (int i=0;i<(int)To_UpdF.size();i++)
    {
        int FIndex=To_UpdF[i];
        for (int j=0;j<(int)Faces[FIndex].size();j++)
        {
            int IndexV=Faces[FIndex][j];
            To_UpdVSet.insert(IndexV);
        }
    }

    std::vector<int> To_UpdVVect(To_UpdVSet.begin(),To_UpdVSet.end());

    for (int i=0;i<(int)To_UpdVVect.size();i++)
    {
        //get the faces to update
        int VIndex=To_UpdVVect[i];
        assert(VIndex>=0);
        assert(VIndex<(int)VF.size());
        for (int j=0;j<(int)To_UpdF.size();j++)
        {
            int FIndex=To_UpdF[j];
            if (!HasVFAdj(VF,VIndex,FIndex))continue;
            VFDetach(VF,VFIndex,VIndex,FIndex);
        }
    }

    //then re-detach
    for (int i=0;i<(int)To_UpdF.size();i++)
    {
        //get the faces to update
        int FIndex=To_UpdF[i];
        assert(FIndex>=0);
        assert(FIndex<(int)Faces.size());
        for (int j=0;j<(int)Faces[FIndex].size();j++)
        {
            int VIndex=Faces[FIndex][j];
            VFAttach(Faces,VF,VFIndex,VIndex,FIndex);
        }
    }
}

void VFDetach(const std::vector<std::vector<int> > &Faces,
              const std::vector<int> &To_UpdF,
              std::vector<std::vector<int> > &VF,
              std::vector<std::vector<int> > &VFIndex)
{
    //first detach all old ones
    for (int i=0;i<(int)To_UpdF.size();i++)
    {
        //get the faces to update
        int FIndex=To_UpdF[i];
        assert(FIndex>=0);
        assert(FIndex<(int)Faces.size());
        for (int j=0;j<(int)Faces[FIndex].size();j++)
        {
            int VIndex=Faces[FIndex][j];
            VFDetach(VF,VFIndex,VIndex,FIndex);
        }
    }
}

void VFAttach(const std::vector<std::vector<int> > &Faces,
              const std::vector<int> &To_UpdF,
              std::vector<std::vector<int> > &VF,
              std::vector<std::vector<int> > &VFIndex)
{
    for (int i=0;i<(int)To_UpdF.size();i++)
    {
        //get the faces to update
        int FIndex=To_UpdF[i];
        assert(FIndex>=0);
        assert(FIndex<(int)Faces.size());
        for (int j=0;j<(int)Faces[FIndex].size();j++)
        {
            int VIndex=Faces[FIndex][j];
            VFAttach(Faces,VF,VFIndex,VIndex,FIndex);
        }
    }
}

//void FFDetach(std::vector<std::vector<int> > &VF,
//              std::vector<std::vector<int> > &VFIndex,
//              const int &VIndex,
//              const int &FIndex)
//{
//}

//void FFAttach(const std::vector<std::vector<int> > &Faces,
//              std::vector<std::vector<int> > &VF,
//              std::vector<std::vector<int> > &VFIndex,
//              const int &VIndex,
//              const int &FIndex)
//{

//}

void UpdateFF(const std::vector<std::vector<int> > &Faces,
              const std::vector<int> &To_UpdF,
              std::vector<std::vector<int> > &NextF,
              std::vector<std::vector<int> > &NextE)
{
    //get the involved edges
    typedef typename std::pair<int,int> intPair;
    std::set<intPair> ToUpdEdges;
    for (size_t i=0;i<To_UpdF.size();i++)
    {
        int IndexF=To_UpdF[i];
        int SizeF=(int)Faces[IndexF].size();
        for (size_t j=0;j<Faces[IndexF].size();j++)
        {
            int IndexV0=Faces[IndexF][j];
            int IndexV1=Faces[IndexF][(j+1)%SizeF];
            intPair KeyE(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            ToUpdEdges.insert(KeyE);
        }
    }

    //get also Old neighbours
    std::set<int> Kernel_F;
    for (size_t i=0;i<To_UpdF.size();i++)
    {
        int IndexF=To_UpdF[i];
        Kernel_F.insert(IndexF);
        for (size_t j=0;j<NextF[IndexF].size();j++)
        {
            if (NextF[IndexF][j]==-1)continue;
            Kernel_F.insert(NextF[IndexF][j]);
        }
    }
    std::vector<int> Kernel_FV(Kernel_F.begin(),Kernel_F.end());

    //get edges to faces, edge I
    std::map<intPair, std::vector<intPair> > EdgeVToFaceE;

    for (size_t i=0;i<Kernel_F.size();i++)
    {
        int IndexF=Kernel_FV[i];
        int SizeF=(int)Faces[IndexF].size();
        for (size_t j=0;j<Faces[IndexF].size();j++)
        {
            int IndexV0=Faces[IndexF][j];
            int IndexV1=Faces[IndexF][(j+1)%SizeF];
            intPair KeyE(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            //check if the one to be updated or not
            if (ToUpdEdges.count(KeyE)==0)continue;

            intPair FaceE(IndexF,j);
            EdgeVToFaceE[KeyE].push_back(FaceE);
        }
    }

    for (size_t i=0;i<Kernel_F.size();i++)
    {
        int IndexF=Kernel_FV[i];
        int SizeF=(int)Faces[IndexF].size();
        for (size_t j=0;j<Faces[IndexF].size();j++)
        {
            int IndexV0=Faces[IndexF][j];
            int IndexV1=Faces[IndexF][(j+1)%SizeF];
            intPair KeyE(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
            if(EdgeVToFaceE.count(KeyE)==0)continue;

            if (EdgeVToFaceE[KeyE].size()==1)
            {
                //in this case is a border
                NextF[IndexF][j]=-1;
                NextE[IndexF][j]=-1;
                continue;
            }
            if (EdgeVToFaceE[KeyE].size()==2)
            {
                int IndexF0=EdgeVToFaceE[KeyE][0].first;
                int IndexE0=EdgeVToFaceE[KeyE][0].second;
                int IndexF1=EdgeVToFaceE[KeyE][1].first;
                int IndexE1=EdgeVToFaceE[KeyE][1].second;

                //in this case is a border
                NextF[IndexF0][IndexE0]=IndexF1;
                NextE[IndexF0][IndexE0]=IndexE1;
                NextF[IndexF1][IndexE1]=IndexF0;
                NextE[IndexF1][IndexE1]=IndexE0;
                continue;
            }
            std::cout<<"SIZE E:"<<EdgeVToFaceE[KeyE].size()<<std::endl;
            //can't be
            assert(0);
        }
    }

}

bool CheckVFAdjacency(const std::vector<std::vector<int> > &Faces,
                      const std::vector<std::vector<int> > &VF,
                      const std::vector<std::vector<int> > &VFIndex)
{
    //check faces to entry
    for (int i=0;i<(int)Faces.size();i++)
        for (int j=0;j<(int)Faces[i].size();j++)
        {
            int CurrVert=Faces[i][j];

            assert(CurrVert>=0);
            assert(CurrVert<(int)VF.size());
            assert(CurrVert<(int)VFIndex.size());

            std::vector<int>::const_iterator IteF=std::find(VF[CurrVert].begin(),VF[CurrVert].end(), i);
            if (IteF==VF[CurrVert].end())return false;
            size_t indexE = std::distance(VF[CurrVert].begin(), IteF);

            assert(indexE>=0);
            assert(indexE<VF[CurrVert].size());
            assert(VF[CurrVert][indexE]==i);

            assert(indexE<VFIndex[CurrVert].size());
            if(VFIndex[CurrVert][indexE]!=j)return false;
        }

    //check entries to faces
    for (int i=0;i<(int)VF.size();i++)
        for (int j=0;j<(int)VF[i].size();j++)
        {
            assert(VF[i].size()==VFIndex[i].size());

            int IndexF=VF[i][j];
            int IndexLocV=VFIndex[i][j];
            assert(IndexF<(int)Faces.size());
            assert(IndexLocV<(int)Faces[IndexF].size());
            if (Faces[IndexF][IndexLocV]!=i)return false;
        }
    return true;
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

int NumSharedFaces(const int &VIndex0,const int &VIndex1,
                   const std::vector<std::vector<int> > &VF)
{
    assert(VIndex0>=0);
    assert(VIndex1>=0);
    assert(VIndex0<(int)VF.size());
    assert(VIndex1<(int)VF.size());
    int SharedE=0;
    std::set<int> SetF0;
    for (size_t i=0;i<VF[VIndex0].size();i++)
        SetF0.insert(VF[VIndex0][i]);

     for (size_t i=0;i<VF[VIndex1].size();i++)
     {
         int IndexF1=VF[VIndex1][i];
         if (SetF0.count(IndexF1)>0)SharedE++;
     }
     return SharedE;
}

}

#endif
