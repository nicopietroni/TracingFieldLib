#ifndef INITIALIZE_GRAPH
#define INITIALIZE_GRAPH

#include <cmath>
#include <vector>
#include <stdio.h>
#include "../point3.h"
#include "../frame_field.h"
#include "../cross_field.h"
#include "../polyhedra_mesh.h"
#include "graph.h"
#include "graph_frame_field.h"
#include "graph_cross_field.h"
#include "connectivity.h"
#include "best_rotation.h"

#define MAX_ANGLE_FOR_DIRECT 90

namespace Geo
{
template <class ScalarType>
class TragingGraphFunctions
{

    inline static ScalarType DriftDistance(const ScalarType &length,const ScalarType &angle,
                                           const ScalarType &max_angle,const ScalarType &drift_penalty)
    {
        ScalarType drift=pow(angle/max_angle,2);
        ScalarType val=length* (ScalarType(1) + drift_penalty * drift);
        return (val);//*Weight;
    }

    static ScalarType DriftAngle(const Geo::Point3<ScalarType> &P0,
                                 const Geo::Point3<ScalarType> &P1,
                                 const Geo::Point3<ScalarType> &D0,
                                 const Geo::Point3<ScalarType> &D1)
    {
        assert(Dot<ScalarType>(D0,D1)>=0);
        Geo::Point3<ScalarType> AvgP=(P0+P1)/(ScalarType)2;
        Geo::Point3<ScalarType> AvgD=(D0+D1)/(ScalarType)2;

        Geo::Point3<ScalarType> DirTest=P1-AvgP;
        ScalarType TestAngle=AngleDeg(DirTest,AvgD);
        return (TestAngle);
    }

    static bool GetVisibilityInfo(const Geo::Point3<ScalarType> &P0,
                                  const Geo::Point3<ScalarType> &P1,
                                  const Geo::Point3<ScalarType> &D0,
                                  const Geo::Point3<ScalarType> &D1,
                                  const ScalarType &drift_penalty,
                                  const ScalarType maxVisAngle,
                                  bool &Visible,
                                  ScalarType &EdgeW)
    {
        Geo::Point3<ScalarType> Dir0=D0;
        Geo::Point3<ScalarType> Dir1=D1;
        Dir0.Normalize();
        Dir1.Normalize();
        ScalarType AngleDir=DriftAngle(P0,P1,Dir0,Dir1);
        assert(maxVisAngle<=MAX_ANGLE_FOR_DIRECT);

        EdgeW=std::numeric_limits<ScalarType>::max();
        Visible=false;

        if (AngleDir>MAX_ANGLE_FOR_DIRECT)
            return false;

        if (AngleDir<=maxVisAngle)
            Visible=true;

        ScalarType Len=(P1-P0).Norm();
        EdgeW=DriftDistance(Len,AngleDir,maxVisAngle,drift_penalty);
        return true;
    }

    static bool GetVisibilityInfoCross(const Field::CrossF<ScalarType> &F0,
                                       const Field::CrossF<ScalarType> &F1,
                                       const int Dir0,const int Dir1,
                                       const ScalarType &drift_penalty,
                                       const ScalarType maxVisAngle,
                                       bool &Visible,
                                       ScalarType &EdgeW)
    {
        assert(Dir0>=0);
        assert(Dir1>=0);
        assert(Dir0<4);
        assert(Dir1<4);

        Geo::Point3<ScalarType> P0=F0.Pos;
        Geo::Point3<ScalarType> P1=F1.Pos;
        Geo::Point3<ScalarType> D0=F0.GetDirection(Dir0);
        Geo::Point3<ScalarType> D1=F1.GetDirection(Dir1);

        return (GetVisibilityInfo(P0,P1,D0,D1,drift_penalty,maxVisAngle,Visible,EdgeW));
    }

    static bool GetVisibilityInfoFrame(const Field::FrameField<ScalarType> &F0,
                                       const Field::FrameField<ScalarType> &F1,
                                       const int Dir0,const int Dir1,
                                       const ScalarType &drift_penalty,
                                       const ScalarType maxVisAngle,
                                       bool &Visible,
                                       ScalarType &EdgeW)
    {
        assert(Dir0>=0);
        assert(Dir1>=0);
        assert(Dir0<6);
        assert(Dir1<6);

        Geo::Point3<ScalarType> P0=F0.Pos;
        Geo::Point3<ScalarType> P1=F1.Pos;
        Geo::Point3<ScalarType> D0=F0.GetDirection(Dir0);
        Geo::Point3<ScalarType> D1=F1.GetDirection(Dir1);

        return (GetVisibilityInfo(P0,P1,D0,D1,drift_penalty,maxVisAngle,Visible,EdgeW));
    }

    static int AddAdjacencyFrame(TracingGraph<ScalarType> &TGraph,
                                 const int IndexN0,
                                 const std::vector<std::vector<int> > &Adjacency,
                                 const std::set<std::pair<int,int> > &DirectConn,
                                 const std::vector<Field::FrameField<ScalarType> > &FrameField,
                                 const ScalarType &drift_penalty,
                                 const ScalarType maxAngle)
    {
        //get the index of polihedron
        int IndexPoly0=Field::FrameFieldGraph::PolyOfNode(IndexN0);
        //and the current direction
        int IndexDir0=Field::FrameFieldGraph::FrameOfNode(IndexN0);

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
            int IndexN1=Field::FrameFieldGraph::NodeOfPolyFrame(IndexPoly1,IndexDir1);


            //CHECK VISIBILITY
            ScalarType EdgeW;
            bool IsVisibleMaxA;
            bool VisibleAbs=GetVisibilityInfoFrame(FrameField[IndexPoly0],
                                                   FrameField[IndexPoly1],
                                                   IndexDir0,IndexDir1,
                                                   drift_penalty,maxAngle,
                                                   IsVisibleMaxA,EdgeW);
            if (!VisibleAbs)
                continue;

            assert(IndexN0!=IndexN1);
            std::pair<int,int> Entry(std::min(IndexPoly0,IndexPoly1),
                                     std::max(IndexPoly0,IndexPoly1));
            bool IsDirect=(DirectConn.count(Entry)>0);
            TGraph.AddNeighbour(IndexN0,IndexN1,EdgeW,
                                true,IsDirect,IsVisibleMaxA);
        }
        return (TGraph.NeighSize(IndexN0));
    }

    static int AddAdjacencyCross(TracingGraph<ScalarType> &TGraph,
                                 const int IndexN0,
                                 const std::vector<std::vector<int> > &Adjacency,
                                 const std::set<std::pair<int,int> > &DirectConn,
                                 const std::vector<Field::CrossF<ScalarType> > &CrossField,
                                 const ScalarType &drift_penalty,
                                 const ScalarType maxAngle)
    {
        //get the index of polihedron
        int IndexVert0=Field::CrossFieldGraph::ElemOfNode(IndexN0);
        //and the current direction
        int IndexDir0=Field::CrossFieldGraph::DirOfNode(IndexN0);

        assert(IndexVert0<(int)Adjacency.size());
        assert(IndexVert0>=0);
        assert(IndexDir0<4);
        assert(IndexDir0>=0);

        //then get all possible neighbours
        for (size_t j=0;j<Adjacency[IndexVert0].size();j++)
        {
            //get the index of polihedron
            int IndexVert1=Adjacency[IndexVert0][j];
            assert(IndexVert1<(int)CrossField.size());
            assert(IndexVert1>=0);

            //and the current direction
            int IndexDir1=CrossField[IndexVert1].FollowDirection(CrossField[IndexVert0],IndexDir0);

            assert(IndexVert1<(int)Adjacency.size());
            assert(IndexDir1>=0);
            assert(IndexDir1<4);
            assert(IndexDir1>=0);

            //then get the new node
            int IndexN1=Field::CrossFieldGraph::NodeOfElemDir(IndexVert1,IndexDir1);
;
            //CHECK VISIBILITY
            ScalarType EdgeW;
            bool IsVisibleMaxA;
            bool VisibleAbs=GetVisibilityInfoCross(CrossField[IndexVert0],
                                                   CrossField[IndexVert1],
                                                   IndexDir0,IndexDir1,
                                                   drift_penalty,maxAngle,
                                                   IsVisibleMaxA,EdgeW);
            if (!VisibleAbs)
                continue;

            assert(IndexN0!=IndexN1);
            std::pair<int,int> Entry(std::min(IndexVert0,IndexVert1),
                                     std::max(IndexVert0,IndexVert1));
            bool IsDirect=(DirectConn.count(Entry)>0);
            TGraph.AddNeighbour(IndexN0,IndexN1,EdgeW,
                                true,IsDirect,IsVisibleMaxA);
        }
        return (TGraph.NeighSize(IndexN0));
        //return (TGraph.N[IndexN0].size());
    }

    static void SetNodeTypeFrame(TracingGraph<ScalarType> &TGraph,
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
            int NodeN=Field::FrameFieldGraph::NodeOfPolyFrame(IndexP,DirN);
            int NodeT0=Field::FrameFieldGraph::NodeOfPolyFrame(IndexP,DirT0);
            int NodeT1=Field::FrameFieldGraph::NodeOfPolyFrame(IndexP,DirT1);
            int NodeT0Opp=Field::FrameFieldGraph::NodeOfPolyFrame(IndexP,DirOppT0);
            int NodeT1Opp=Field::FrameFieldGraph::NodeOfPolyFrame(IndexP,DirOppT1);

            TGraph.NodeType(NodeN)=ExitDomain;
            TGraph.NodeType(NodeT0)=TangentDomain;
            TGraph.NodeType(NodeT1)=TangentDomain;
            TGraph.NodeType(NodeT0Opp)=TangentDomain;
            TGraph.NodeType(NodeT1Opp)=TangentDomain;
        }

    }

    static void GetVertParallelExitDir(const std::vector<std::vector<int> > &FaceConnectivity,
                                       const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                       std::vector<std::vector<Geo::Point3<ScalarType> > > &BorderVDirs,
                                       std::vector<Geo::Point3<ScalarType> > &BorderExitVDirs)
    {
        ScalarType eps=0.0000001;
        //update Normals
        std::vector<Point3<ScalarType> > FaceNormals;
        std::vector<Point3<ScalarType> > VertNormals;
        Geo::ComputeNormals(VertPos,FaceConnectivity,FaceNormals,VertNormals);

        std::vector<std::vector<int> > NextF,NextE;
        Geo::getFFAdjacency(FaceConnectivity,NextF,NextE);

        BorderVDirs.resize(VertPos.size());
        BorderExitVDirs.resize(VertPos.size(),Geo::Point3<ScalarType>(0,0,0));

        for (size_t i=0;i<NextF.size();i++)
            for (size_t j=0;j<NextF[i].size();j++)
            {
                if (NextF[i][j]!=-1)continue;//only border
                int sizeF=FaceConnectivity[i].size();
                int IndexV0=FaceConnectivity[i][j];
                int IndexV1=FaceConnectivity[i][(j+1)%sizeF];
                assert(IndexV0!=IndexV1);
                assert(IndexV0>=0);
                assert(IndexV1>=0);
                assert(IndexV0<(int)VertPos.size());
                assert(IndexV1<(int)VertPos.size());
                Geo::Point3<ScalarType> Pos0=VertPos[IndexV0];
                Geo::Point3<ScalarType> Pos1=VertPos[IndexV1];
                Geo::Point3<ScalarType> Dir=Pos1-Pos0;
                Dir.Normalize();

                //bring the vector in the right space
                Geo::Point3<ScalarType> NormF=FaceNormals[i];
                Geo::Point3<ScalarType> Norm0=VertNormals[IndexV0];
                Geo::Point3<ScalarType> Norm1=VertNormals[IndexV1];
                assert(fabs(NormF.Norm()-1)<eps);
                assert(fabs(Norm0.Norm()-1)<eps);
                assert(fabs(Norm1.Norm()-1)<eps);
                Eigen::Matrix3<ScalarType> Rot0=Geo::PointMatch<ScalarType>::BestRotationMatrix(NormF,Norm0);
                Eigen::Matrix3<ScalarType> Rot1=Geo::PointMatch<ScalarType>::BestRotationMatrix(NormF,Norm0);

                //the direction should point on the other side in this case for IndexV1
                Geo::Point3<ScalarType> Dir0=Rot0*Dir;
                Geo::Point3<ScalarType> Dir1=Rot1*(Dir);
                BorderVDirs[IndexV0].push_back(Dir0);
                BorderVDirs[IndexV1].push_back(-Dir1);

                //then compute the opposite direction, but should always turn the same side
                Geo::Point3<ScalarType> Opp0=Cross(Dir0,NormF);
                Geo::Point3<ScalarType> Opp1=Cross(Dir1,NormF);
                Opp0.Normalize();
                Opp1.Normalize();
                BorderExitVDirs[IndexV0]+=Opp0;
                BorderExitVDirs[IndexV1]+=Opp1;
            }

        for (size_t i=0;i<BorderExitVDirs.size();i++)
            BorderExitVDirs[i].Normalize();
    }

    static void SetNodeTypeVertCross(TracingGraph<ScalarType> &TGraph,
                                     const std::vector<std::vector<int> > &FaceConnectivity,
                                     const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                     const std::vector<Field::CrossF<ScalarType> > &CrossField)
    {
        std::vector<std::vector<Geo::Point3<ScalarType> > > BorderVDirs;
        std::vector<Geo::Point3<ScalarType> > BorderExitVDirs;
        GetVertParallelExitDir(FaceConnectivity,VertPos,BorderVDirs,BorderExitVDirs);

        assert(BorderVDirs.size()==VertPos.size());
        assert(BorderExitVDirs.size()==VertPos.size());
        assert(CrossField.size()==VertPos.size());
        for (size_t i=0;i<VertPos.size();i++)
        {
            if (BorderVDirs[i].size()==0)continue;

            //set the exit dir
            assert(BorderExitVDirs[i]!=Geo::Point3<ScalarType>(0,0,0));
            int ExitDirIdx=CrossField[i].BestDirectionI(BorderExitVDirs[i]);
            assert(ExitDirIdx>=0);
            assert(ExitDirIdx<4);
            int IndexNOpp=Field::CrossFieldGraph::NodeOfElemDir(i,ExitDirIdx);
            TGraph.NodeType(IndexNOpp)=ExitDomain;

            for (size_t j=0;j<BorderVDirs[i].size();j++)
            {
                assert(BorderVDirs[i][j]!=Geo::Point3<ScalarType>(0,0,0));
                int BorderDirIdx=CrossField[i].BestDirectionI(BorderVDirs[i][j]);
                assert(BorderDirIdx>=0);
                assert(BorderDirIdx<4);
                int IndexNBord=Field::CrossFieldGraph::NodeOfElemDir(i,BorderDirIdx);
                TGraph.NodeType(IndexNBord)=TangentDomain;
            }
        }
    }

    static void FixNonExitSink(TracingGraph<ScalarType> &TGraph,
                               const std::vector<std::vector<int> > &Adjacency,
                               const std::set<std::pair<int,int> > &EdgeSet,
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
                AddedN=AddAdjacencyFrame(TGraph,(int)i,Adjacency,EdgeSet,FrameField,drift_penalty,CurrA);
            }while ((AddedN==0)&&(CurrA<270));
            if (AddedN==0)
                std::cout<<"WARNING INTERNAL SINK"<<std::endl;
        }
    }

    static void RemoveBorderExitNeigh(TracingGraph<ScalarType> &TGraph,
                                      bool RemoveTangent=true)
    {
        for (size_t i=0;i<TGraph.NumNodes();i++)
        {
            if (TGraph.NodeType(i)==ExitDomain)
                TGraph.RemoveNeighbours(i);
            if ((RemoveTangent)&&(TGraph.NodeType(i)==TangentDomain))
                TGraph.RemoveNeighbours(i);
        }
    }

    static void InitAdjacencyFromPolyFrameField(TracingGraph<ScalarType> &TGraph,
                                                const std::vector<std::vector<int> > &Adjacency,
                                                const std::set<std::pair<int,int> > &EdgeSet,
                                                const std::vector<Field::FrameField<ScalarType> > &FrameField,
                                                const ScalarType &drift_penalty,
                                                const ScalarType maxAngle)
    {

        TGraph.AllocateNodes(FrameField.size()*6);
        //for each node
        for (size_t i=0;i<TGraph.NumNodes();i++)
            AddAdjacencyFrame(TGraph,(int)i,Adjacency,EdgeSet,FrameField,drift_penalty,maxAngle);

        TGraph.SortAllNeigh();
    }

    static void InitAdjacencyFromVertCrossField(TracingGraph<ScalarType> &TGraph,
                                                const std::vector<std::vector<int> > &Adjacency,
                                                const std::set<std::pair<int,int> > &EdgeSet,
                                                const std::vector<Field::CrossF<ScalarType> > &CrossField,
                                                const ScalarType &drift_penalty,
                                                const ScalarType maxAngle)
    {

        //instanciate 6 nodes per each frame (one for each direction)
        TGraph.AllocateNodes(CrossField.size()*4);

        //for each node
        for (size_t i=0;i<TGraph.NumNodes();i++)
            AddAdjacencyCross(TGraph,(int)i,Adjacency,EdgeSet,CrossField,drift_penalty,maxAngle);

        TGraph.SortAllNeigh();
    }

    static void GetEdgeSetFromAdj(const std::vector<std::vector<int> > &Adjacency,
                                  std::set<std::pair<int,int> > &EdgeSet)
    {
        EdgeSet.clear();
        for (size_t i=0;i<Adjacency.size();i++)
            for (size_t j=0;j<Adjacency[i].size();j++)
            {
                int IndexN0=i;
                int IndexN1=Adjacency[i][j];
                std::pair<int,int> Key(std::min(IndexN0,IndexN1),
                                       std::max(IndexN0,IndexN1));
                EdgeSet.insert(Key);
            }
    }

    static void GetCrossAdjDirectSet(const std::vector<std::vector<int> > &Faces,
                                     std::set<std::pair<int,int> > &EdgeSet)
    {
        std::vector<std::vector<int> > Adjacency0;
        Geo::FindVertexStar(Faces,Adjacency0,0);
        GetEdgeSetFromAdj(Adjacency0,EdgeSet);
    }

    static void GetFrameAdjDirectSet(const std::vector<std::vector<int> > &PolyConnectivity,
                                     std::set<std::pair<int,int> > &EdgeSet)
    {
        std::vector<std::vector<int> > Adjacency0;
        Tetra::GetPPKernelAlongVertices(PolyConnectivity,Adjacency0,0);
        GetEdgeSetFromAdj(Adjacency0,EdgeSet);
    }

    static bool IsProperConnectivity(const std::vector<std::vector<int> > &Adjacency)
    {
        //no link to itself
        for (size_t i=0;i<Adjacency.size();i++)
            for (size_t j=0;j<Adjacency[i].size();j++)
               if (Adjacency[i][j]==i)
               {std::cout<<"DE0"<<std::endl;return false;}

        //no double
        for (size_t i=0;i<Adjacency.size();i++)
            for (size_t j=0;j<Adjacency[i].size()-1;j++)
                 for (size_t k=(j+1);k<Adjacency[i].size();k++)
                    if (Adjacency[i][j]==Adjacency[i][k])
                    {std::cout<<"DE1"<<std::endl;return false;}

        return true;
    }

public:

    static void InitFromVertCrossField(TracingGraph<ScalarType> &TGraph,
                                       const std::vector<Geo::Point3<ScalarType> > &VertPos,
                                       const std::vector<std::vector<int> > &Faces,
                                       const std::vector<Field::CrossF<ScalarType> > &CrossField,
                                       const ScalarType &drift_penalty=100,
                                       const ScalarType maxAngle=45,
                                       const size_t diffuse_step=0)
    {
        TGraph.Clear();
        assert(VertPos.size()==CrossField.size());
        std::vector<std::vector<int> > Adjacency;
        Geo::FindVertexStar(Faces,Adjacency,diffuse_step);

        std::set<std::pair<int,int> > EdgeSet;
        GetCrossAdjDirectSet(Faces,EdgeSet);

        InitAdjacencyFromVertCrossField(TGraph,Adjacency,EdgeSet,CrossField,drift_penalty,maxAngle);

        SetNodeTypeVertCross(TGraph,Faces,VertPos,CrossField);
        RemoveBorderExitNeigh(TGraph);


        TGraph.InitTracing();
        TGraph.SortAllNeigh();
    }

    static void InitFromPolyFrameField(TracingGraph<ScalarType> &TGraph,
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
        //GetGraphComplexAdjacency(VertPos,PolyConnectivity,Adjacency,diffuse_step);
        Tetra::GetPPKernelAlongVertices(VertPos.size(),PolyConnectivity,Adjacency,diffuse_step);

        assert(IsProperConnectivity(Adjacency));

        std::set<std::pair<int,int> > EdgeSet;
        GetCrossAdjDirectSet(PolyConnectivity,EdgeSet);

        InitAdjacencyFromPolyFrameField(TGraph,Adjacency,EdgeSet,FrameField,drift_penalty,maxAngle);

        SetNodeTypeFrame(TGraph,PolyConnectivity,VertPos,FrameField);

        //RemoveBorderPolyField(TGraph);
        RemoveBorderExitNeigh(TGraph,false);

        assert(TGraph.ProperNeigh());

        //FixNonExitSink(TGraph,Adjacency,EdgeSet,FrameField,drift_penalty,maxAngle);

        TGraph.InitTracing();

        TGraph.SortAllNeigh();

        assert(TGraph.ProperNeigh());
    }
};
}
#endif
