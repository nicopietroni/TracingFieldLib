#ifndef INITIALIZE_GRAPH
#define INITIALIZE_GRAPH

#include <cmath>
#include <vector>
#include <stdio.h>
#include "../point3.h"
#include "../frame_field.h"
#include "../polyhedra_mesh.h"
#include "graph.h"
#include "graph_frame_field.h"

namespace Geo
{
template <class ScalarType>
class TragingGraphFunctions
{


    inline static ScalarType DriftDistance(const ScalarType &length,const ScalarType &angle,
                                           const ScalarType &max_angle,const ScalarType &drift_penalty)
    {
        //    if (max_angle == ScalarType(0))
        //        return (length*Weight);
        //return length * (ScalarType(1) + drift_penalty * angle / max_angle)*Weight;
        ScalarType drift=pow(angle/max_angle,2);
        ScalarType val=length* (ScalarType(1) + drift_penalty * drift);
        return (val);//*Weight;
    }

    static bool IsVisible(const Field::FrameField<ScalarType> &F0,
                          const Field::FrameField<ScalarType> &F1,
                          const int Dir0,const int Dir1,
                          const ScalarType &drift_penalty,
                          const ScalarType maxAngle,
                          ScalarType &EdgeW)
    {
        assert(Dir0>=0);
        assert(Dir1>=0);
        assert(Dir0<6);
        assert(Dir1<6);

        EdgeW=-1;
        Geo::Point3<ScalarType> P0=F0.Pos;
        Geo::Point3<ScalarType> P1=F1.Pos;
        Geo::Point3<ScalarType> D0=F0.GetDirection(Dir0);
        Geo::Point3<ScalarType> D1=F1.GetDirection(Dir1);
        assert(Dot<ScalarType>(D0,D1)>=0);

        Geo::Point3<ScalarType> AvgP=(P0+P1)/(ScalarType)2;
        Geo::Point3<ScalarType> AvgD=(D0+D1)/(ScalarType)2;

        Geo::Point3<ScalarType> DirTest=P1-AvgP;
        ScalarType TestAngle=AngleDeg(DirTest,AvgD);
        if (TestAngle>maxAngle)
            return false;

        ScalarType Len=(P1-P0).Norm();
        EdgeW=DriftDistance(Len,TestAngle,maxAngle,drift_penalty);
        return true;
    }

    static int AddAdjacency(TracingGraph<ScalarType> &TGraph,
                            const int IndexN0,
                            const std::vector<std::vector<int> > &Adjacency,
                            const std::vector<Field::FrameField<ScalarType> > &FrameField,
                            const ScalarType &drift_penalty,
                            const ScalarType maxAngle)
    {
        //get the index of polihedron
        int IndexPoly0=PolyOfNode(IndexN0);
        //and the current direction
        int IndexDir0=FrameOfNode(IndexN0);

        assert(IndexPoly0<Adjacency.size());
        assert(IndexPoly0>=0);
        assert(IndexDir0<6);
        assert(IndexDir0>=0);

        //then get all possible neighbours
        for (size_t j=0;j<Adjacency[IndexPoly0].size();j++)
        {
            //get the index of polihedron
            int IndexPoly1=Adjacency[IndexPoly0][j];
            assert(IndexPoly1<FrameField.size());
            assert(IndexPoly1>=0);

            //and the current direction
            int IndexDir1=FrameField[IndexPoly1].FollowDirection(FrameField[IndexPoly0],IndexDir0);

            assert(IndexPoly1<Adjacency.size());
            assert(IndexPoly1>=0);
            assert(IndexDir1<6);
            assert(IndexDir1>=0);

            //then get the new node
            int IndexN1=NodeOfPolyFrame(IndexPoly1,IndexDir1);


            //CHECK VISIBILITY
            ScalarType EdgeW;
            if (!IsVisible(FrameField[IndexPoly0],FrameField[IndexPoly1],
                           IndexDir0,IndexDir1,drift_penalty,maxAngle,EdgeW))
                continue;

            TGraph.N[IndexN0].push_back(IndexN1);
            TGraph.W[IndexN0].push_back(EdgeW);
        }
        return (TGraph.N[IndexN0].size());
    }

    static void SetNodeType(TracingGraph<ScalarType> &TGraph,
                            const std::vector<std::vector<int> > &PolyConnectivity,
                            const std::vector<Geo::Point3<ScalarType> > &VertPos,
                            const std::vector<Field::FrameField<ScalarType> > &FrameField)
    {
        std::vector<Geo::Point3<ScalarType> > VertNormals;

        std::vector<std::vector<int> > ExtFace;
        std::vector<int> ExtPoly;
        std::vector<int> ExtPolyFaceIndex;
        //std::vector<std::vector<std::pair<int,int> > > PPAdjacency;
        std::vector<std::vector<int> > NextP,NextF;
        Tetra::GetPPAdjacency(PolyConnectivity,NextP,NextF);
        Tetra::FFExternalMesh(PolyConnectivity,NextP,NextF,ExtFace,ExtPoly,ExtPolyFaceIndex);

        //    TGraph.Sink.clear();
        //    TGraph.Sink.resize(TGraph.N.size(),false);

        std::vector<Geo::Point3<ScalarType> > FaceNormals;
        ComputeNormals(VertPos,ExtFace,FaceNormals,VertNormals);
        assert(ExtPoly.size()==FaceNormals.size());

        for (size_t i=0;i<ExtPoly.size();i++)
        {
            int IndexP=ExtPoly[i];
            Geo::Point3<ScalarType> Norm=FaceNormals[i];
            assert(IndexP>=0);
            assert(IndexP<FrameField.size());

            int DirN,DirT0,DirT1;
            FrameField[IndexP].GetFrameDir(Norm,DirN,DirT0,DirT1);
            assert(DirN>=0);
            assert(DirN<6);

            int DirOppT0,DirOppT1;
            DirOppT0=Field::OppositeFrameDir(DirT0);
            DirOppT1=Field::OppositeFrameDir(DirT1);

            //then remove negibours
            int NodeN=NodeOfPolyFrame(IndexP,DirN);
            int NodeT0=NodeOfPolyFrame(IndexP,DirT0);
            int NodeT1=NodeOfPolyFrame(IndexP,DirT1);
            int NodeT0Opp=NodeOfPolyFrame(IndexP,DirOppT0);
            int NodeT1Opp=NodeOfPolyFrame(IndexP,DirOppT1);

            TGraph.NodeType(NodeN)=ExitDomain;
            TGraph.NodeType(NodeT0)=TangentDomain;
            TGraph.NodeType(NodeT1)=TangentDomain;
            TGraph.NodeType(NodeT0Opp)=TangentDomain;
            TGraph.NodeType(NodeT1Opp)=TangentDomain;
        }
    }

    static void FixNonExitSink(TracingGraph<ScalarType> &TGraph,
                               const std::vector<std::vector<int> > &Adjacency,
                               const std::vector<Field::FrameField<ScalarType> > &FrameField,
                               const ScalarType &drift_penalty,
                               const ScalarType maxAngle)
    {
        for (size_t i=0;i<TGraph.NumNodes();i++)
        {
            if (TGraph.NodeType(i)==ExitDomain)continue;
            if (!TGraph.IsSink(i))continue;
            int AddedN=0;
            ScalarType CurrA=maxAngle;
            do{
                CurrA+=10;
                AddedN=AddAdjacency(TGraph,(int)i,Adjacency,FrameField,drift_penalty,CurrA);
            }while ((AddedN==0)&&(CurrA<270));
            if (AddedN==0)
                std::cout<<"WARNING INTERNAL SINK"<<std::endl;
        }
    }

    static void RemoveBorderPolyField(TracingGraph<ScalarType> &TGraph)
    {
        for (size_t i=0;i<TGraph.NumNodes();i++)
        {
            if (TGraph.NodeType(i)==ExitDomain)//||(TGraph.NodeType(i)==TangentDomain))
                TGraph.RemoveNeighbours(i);
        }
    }

    static void InitAdjacencyFromPolyField(TracingGraph<ScalarType> &TGraph,
                                           const std::vector<std::vector<int> > &Adjacency,
                                           const std::vector<Field::FrameField<ScalarType> > &FrameField,
                                           const ScalarType &drift_penalty,
                                           const ScalarType maxAngle)
    {

        //instanciate 6 nodes per each frame (one for each direction)
        TGraph.N.resize(FrameField.size()*6);
        TGraph.W.resize(FrameField.size()*6);
        TGraph.NT.clear();
        TGraph.NT.resize(TGraph.N.size(),Internal);

        //for each node
        for (size_t i=0;i<TGraph.N.size();i++)
            AddAdjacency(TGraph,(int)i,Adjacency,FrameField,drift_penalty,maxAngle);

        TGraph.SortAllNeighbyW();
    }

public:

    static void GetGraphComplexAdjacency(const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                         const std::vector<std::vector<int> > &PolyConnectivity,
                                         std::vector<std::vector<int> > &Adjacency,
                                         size_t diffuse_step=0)
    {
        assert(diffuse_step==0);//to be implemented

        std::vector<std::vector<std::pair<int,int> > > VPAdjacency;
        Tetra::GetVPAdjacency(VertPos.size(),PolyConnectivity,VPAdjacency);

        Adjacency.clear();
        Adjacency.resize(PolyConnectivity.size());

        for (size_t i=0;i<VPAdjacency.size();i++)
        {
            std::vector<int>  Polys;
            for (size_t j=0;j<VPAdjacency[i].size();j++)
                Polys.push_back(VPAdjacency[i][j].first);

            for (size_t k=0;k<Polys.size();k++)
                for (size_t w=0;w<Polys.size();w++)
                {
                    if (k==w)continue;
                    int IndexP0=Polys[k];
                    int IndexP1=Polys[w];
                    assert(IndexP0!=IndexP1);
                    assert(IndexP0>=0);
                    assert(IndexP1>=0);
                    assert(IndexP0<Adjacency.size());
                    assert(IndexP1<Adjacency.size());
                    Adjacency[IndexP0].push_back(IndexP1);
                    Adjacency[IndexP1].push_back(IndexP0);
                }
        }

        //remove duplicates
        for (size_t i=0;i<Adjacency.size();i++)
        {
            std::sort(Adjacency[i].begin(),Adjacency[i].end());
            std::vector<int>::iterator last = std::unique(Adjacency[i].begin(),Adjacency[i].end());
            Adjacency[i].erase(last, Adjacency[i].end());
        }
    }

    static void InitFromPolyField(TracingGraph<ScalarType> &TGraph,
                                  const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                  const std::vector<std::vector<int> > &PolyConnectivity,
                                  const std::vector<Field::FrameField<ScalarType> > &FrameField,
                                  const ScalarType &drift_penalty=100,
                                  const ScalarType maxAngle=45,
                                  const size_t diffuse_step=0)
    {
        TGraph.Clear();

        assert(PolyConnectivity.size()==FrameField.size());
        std::vector<std::vector<int> > Adjacency;
        GetGraphComplexAdjacency(VertPos,PolyConnectivity,Adjacency,diffuse_step);
        InitAdjacencyFromPolyField(TGraph,Adjacency,FrameField,drift_penalty,maxAngle);

        SetNodeType(TGraph,PolyConnectivity,VertPos,FrameField);
        RemoveBorderPolyField(TGraph);

        FixNonExitSink(TGraph,Adjacency,FrameField,drift_penalty,maxAngle);

        TGraph.InitTracing();
        TGraph.SortAllNeighbyW();
    }
};
}
#endif
