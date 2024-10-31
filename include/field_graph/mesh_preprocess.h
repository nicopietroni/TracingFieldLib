#ifndef TRACING_MESH_PREPROCESS
#define TRACING_MESH_PREPROCESS

#include "field_graph/graph_cross_field.h"
#include "point3.h"
#include <cassert>
#include <connectivity.h>
#include <cstddef>
#include <iostream>
#include <set>
#include <mesh_refine.h>
#include <sys/_types/_size_t.h>
#include <utility>
#include <vector>

namespace Geo
{

    template <class ScalarType>
    class TracingPreprocess
    {
        static void GetAllEdgesOfFaces(const std::vector<std::vector<int>> &Faces,
                                       const std::vector<int> &FacesIdx,
                                       std::vector<std::pair<int, int>> &EdgesV)
        {
            EdgesV.clear();
            for (size_t i = 0; i < FacesIdx.size(); i++)
            {
                int IndexF = FacesIdx[i];
                assert(Faces[IndexF].size() == 3);

                for (size_t j = 0; j < 3; j++)
                {
                    int IndexV0 = Faces[IndexF][j];
                    int IndexV1 = Faces[IndexF][(j + 1) % 3];
                    std::pair<int, int> Entry(std::min(IndexV0, IndexV1), std::max(IndexV0, IndexV1));
                    EdgesV.push_back(Entry);
                }
            }
            std::sort(EdgesV.begin(), EdgesV.end());
            EdgesV.erase(std::unique(EdgesV.begin(), EdgesV.end()), EdgesV.end());
        }

        static void GetAllVertOfFaces(const std::vector<std::vector<int>> &Faces,
                                      const std::vector<int> &FacesIdx,
                                      std::vector<int> &IndexV)
        {
            IndexV.clear();
            for (size_t i = 0; i < FacesIdx.size(); i++)
            {
                int IndexF = FacesIdx[i];
                for (size_t k = 0; k < Faces[IndexF].size(); k++)
                    IndexV.push_back(Faces[IndexF][k]);
            }
            std::sort(IndexV.begin(), IndexV.end());
            IndexV.erase(std::unique(IndexV.begin(), IndexV.end()), IndexV.end());
        }

        static void GetCloseSingFaces(const std::vector<std::vector<int>> &Faces,
                                      const std::vector<int> &SingIndex,
                                      const int extend_step,
                                      std::vector<int> &CloseSingF)
        {
            CloseSingF.clear();

            std::vector<std::vector<int>> VF, VFIndex;
            Geo::getVFAdjacency(Faces, VF, VFIndex);

            std::set<int> SingSet(SingIndex.begin(), SingIndex.end());

            for (size_t i = 0; i < SingIndex.size(); i++)
            {
                int CurrV = SingIndex[i];
                std::vector<int> FacesStar;
                getVFStar(Faces, VF, CurrV, extend_step, FacesStar);

                bool found = false;
                std::vector<int> IndexV;
                GetAllVertOfFaces(Faces, FacesStar, IndexV);
                // int NumSing = 0;
                for (size_t j = 0; j < IndexV.size(); j++)
                {
                    if (IndexV[j] == CurrV)
                        continue;
                    if (SingSet.count(IndexV[j]) == 0)
                        continue;
                    found = true;
                    break;
                }
                if (found)
                    CloseSingF.insert(CloseSingF.end(),
                                      FacesStar.begin(),
                                      FacesStar.end());
            }
        }

        // static bool RefineCloseSingFacesStep(std::vector<Geo::Point3<ScalarType>> &Pos,
        //                                      std::vector<std::vector<int>> &Faces,
        //                                      std::vector<Field::CrossF<ScalarType> > &FaceCurv,
        //                                      std::vector<std::pair<int,int>> &Features,
        //                                      const std::vector<int> &SingIndex,
        //                                      const int extend_step)
        // {
        //     //std::set<std::pair<int,int>> FeatureSet(Features.begin(),Features.end());
        //     std::set<std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>>> MarkedE;
        //     for (size_t i=0;i<Features.size();i++)
        //     {
        //         int IndexV0=Features[i].first;
        //         int IndexV1=Features[i].second;
        //         assert(IndexV0>=0);
        //         assert(IndexV1>=0);
        //         assert(IndexV0<Pos.size());
        //         assert(IndexV1<Pos.size());
        //         Geo::Point3<ScalarType> P0=Pos[IndexV0];
        //         Geo::Point3<ScalarType> P1=Pos[IndexV1];
        //         std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> Key(std::min(P0,P1),std::max(P0,P1));
        //         MarkedE.insert(Key);
        //     }
            
        //     std::vector<int> CloseSingF;
        //     GetCloseSingFaces(Faces, SingIndex, extend_step, CloseSingF);

        //     // then set the edges to be split
        //     std::vector<std::pair<int, int>> To_splitE;
        //     GetAllEdgesOfFaces(Faces, CloseSingF, To_splitE);

        //     if (To_splitE.size() == 0)
        //         return false;

        //     std::vector<ScalarType> Alpha(To_splitE.size(), 0.5);
        //     std::vector<int> OldF;
        //     RefineMesh<ScalarType>::SplitTriangleFaces(Pos, Faces, To_splitE, Alpha,MarkedE,false, OldF);
        //     assert(OldF.size() == Faces.size());

        //     //update cross field
        //     std::vector<Field::CrossF<ScalarType>> OldFaceCurv = FaceCurv;
        //     FaceCurv.clear();
        //     for (size_t i = 0; i < OldF.size(); i++)
        //     {
        //         // int IndexFNew=i;
        //         int IndexFOld = OldF[i];
        //         FaceCurv.push_back(OldFaceCurv[IndexFOld]);

        //         int IndexV0 = Faces[i][0];
        //         int IndexV1 = Faces[i][1];
        //         int IndexV2 = Faces[i][2];
        //         Geo::Point3<ScalarType> P0 = Pos[IndexV0];
        //         Geo::Point3<ScalarType> P1 = Pos[IndexV1];
        //         Geo::Point3<ScalarType> P2 = Pos[IndexV2];
        //         Geo::Point3<ScalarType> BaryF = (P0 + P1 + P2) / 3;
        //         FaceCurv.back().Pos = BaryF;
        //     }

        //     //update features
        //     Features.clear();
        //     for (size_t i = 0; i < Faces.size(); i++)
        //     for (size_t j = 0; j < Faces[i].size(); j++)
        //     {
        //         int IndexV0 = Faces[i][j];
        //         int IndexV1 = Faces[i][(j+1)%Faces[i].size()];
        //         Geo::Point3<ScalarType> P0=Pos[IndexV0];
        //         Geo::Point3<ScalarType> P1=Pos[IndexV1];
        //         std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> Key(std::min(P0,P1),std::max(P0,P1));
        //         //std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> Key1(std::max(P0,P1),std::min(P0,P1));
        //         //if ((MarkedE.count(Key0)==0)&&(MarkedE.count(Key1)==0))continue;
        //         if (MarkedE.count(Key)==0)continue;
        //         Features.push_back(std::pair<int,int>(std::min(IndexV0,IndexV1),
        //                                               std::max(IndexV0,IndexV1)));
        //     }

        //     return true;
        // }

        static bool RefineCloseSingFacesStep(std::vector<Geo::Point3<ScalarType>> &Pos,
                                             std::vector<std::vector<int>> &Faces,
                                             std::vector<Field::CrossF<ScalarType> > &FaceCurv,
                                             std::vector<std::pair<int,int>> &Features,
                                             const std::vector<int> &SingIndex,
                                             const int extend_step)
        {
            
            std::vector<int> CloseSingF;
            GetCloseSingFaces(Faces, SingIndex, extend_step, CloseSingF);

            // then set the edges to be split
            std::vector<std::pair<int, int>> To_splitE;
            GetAllEdgesOfFaces(Faces, CloseSingF, To_splitE);

            if (To_splitE.size() == 0)
                return false;

            std::vector<Geo::Point3<ScalarType>> NewPos;
            for (size_t i=0;i<To_splitE.size();i++)
            {
                int IndexV0=To_splitE[i].first;
                int IndexV1=To_splitE[i].second;
                Geo::Point3<ScalarType> P0=Pos[IndexV0];
                Geo::Point3<ScalarType> P1=Pos[IndexV1];
                if (P1<P0)std::swap(P0,P1);
                Geo::Point3<ScalarType> MidP=P0*0.5+P1*0.5;
                NewPos.push_back(MidP);
            }
            std::vector<int> OldF;
            //RefineMesh<ScalarType>::SplitTriangleFaces(Pos, Faces, To_splitE, Alpha,MarkedE,false, OldF);
            RefineMesh<ScalarType>::SplitTriangleFaces(Pos,Faces,To_splitE,NewPos,Features,OldF);
            assert(OldF.size() == Faces.size());

            //update cross field
            std::vector<Field::CrossF<ScalarType>> OldFaceCurv = FaceCurv;
            FaceCurv.clear();
            for (size_t i = 0; i < OldF.size(); i++)
            {
                // int IndexFNew=i;
                int IndexFOld = OldF[i];
                FaceCurv.push_back(OldFaceCurv[IndexFOld]);

                int IndexV0 = Faces[i][0];
                int IndexV1 = Faces[i][1];
                int IndexV2 = Faces[i][2];
                Geo::Point3<ScalarType> P0 = Pos[IndexV0];
                Geo::Point3<ScalarType> P1 = Pos[IndexV1];
                Geo::Point3<ScalarType> P2 = Pos[IndexV2];
                Geo::Point3<ScalarType> BaryF = (P0 + P1 + P2) / 3;
                FaceCurv.back().Pos = BaryF;
            }

            return true;
        }

    public:

        static void RefineCloseSingFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                                         std::vector<std::vector<int>> &Faces,
                                         std::vector<Field::CrossF<ScalarType>> &FaceCurv,
                                         std::vector<std::pair<int,int>> &Features,
                                         const std::vector<int> &SingIndex,
                                         const int extend_step)
        {
            RefineCloseSingFacesStep(Pos, Faces, FaceCurv,Features, SingIndex, extend_step);
            // while (RefineCloseSingFacesStep(Pos, Faces, FaceCurv,Features, SingIndex, extend_step))
            // {}
        }

        static void DesambiguateCoincidentBorders(std::vector<Geo::Point3<ScalarType>> &Pos,
                                             std::vector<std::vector<int>> &Faces,
                                             std::vector<int> *ModVert=NULL,
                                             std::vector<Geo::Point3<ScalarType>> *OrigPos=NULL)
        {
            ScalarType delta= AvgEdgeLen(Pos,Faces)*0.01;
            //ScalarType delta=std::numeric_limits<float>::epsilon()*(ScalarType)2;
            DesambiguateSamePosVertices<ScalarType>(Pos, Faces,delta,ModVert,OrigPos);
        }

        static void SplitFeatures(std::vector<Geo::Point3<ScalarType>> &Pos,
                                  std::vector<std::vector<int>> &Faces,
                                  std::vector<Field::CrossF<ScalarType>> &VertCurv,
                                  std::vector<int> &SingIndex, 
                                  std::vector<int> &SingValue,
                                  std::vector<std::pair<int,int>> &Features)
        {
            
            assert(VertCurv.size()==Pos.size());
            assert(SingIndex.size()==SingValue.size());

            std::map<Geo::Point3<ScalarType>,Field::CrossF<ScalarType> > PosCross;
            for (size_t i=0;i<Pos.size();i++)
                PosCross[Pos[i]]=VertCurv[i];
            
            std::map<Geo::Point3<ScalarType>,int > PosSing;
            for (size_t i=0;i<SingIndex.size();i++)
            {
                Geo::Point3<ScalarType> SingPos=Pos[SingIndex[i]];
                PosSing[SingPos]=SingValue[i];
            }

            std::set<std::pair<Geo::Point3<ScalarType>,Geo::Point3<ScalarType> > > FeaturePos;
            for (size_t i=0;i<Features.size();i++)
            {
               int IndexV0=Features[i].first;
               int IndexV1=Features[i].second;
               Geo::Point3<ScalarType> P0=Pos[IndexV0];
               Geo::Point3<ScalarType> P1=Pos[IndexV1];
               std::pair<Geo::Point3<ScalarType>,Geo::Point3<ScalarType> > KeyF(std::min(P0,P1),std::max(P0,P1));
               FeaturePos.insert(KeyF);
            }

            // std::vector<int> &SingIndex, 
            //                       std::vector<int> &SingValue,
            //std::vector<std::pair<int,int>> SplitE=Features;
            // SplitMeshFromEdges(std::vector<Geo::Point3<ScalarType>> &Pos,
            //             std::vector<std::vector<int>> &Faces,
            //             const std::vector<std::pair<int, int>> &SplitEdges,
            //             std::vector<std::pair<int, int>> &Features)

            SplitMeshFromEdges<ScalarType>(Pos,Faces,Features);//,Features);
            RemoveUnfererencedVert<ScalarType>(Pos,Faces);

            VertCurv.clear();
            for (size_t i=0;i<Pos.size();i++)
            {
                assert(PosCross.count(Pos[i])>0);
                VertCurv.push_back(PosCross[Pos[i]]);
            }
            
            SingIndex.clear();
            SingValue.clear();
            for (size_t i=0;i<Pos.size();i++)
            {
                Geo::Point3<ScalarType> TestPos=Pos[i];
                if (PosSing.count(TestPos)==0)continue;
                SingIndex.push_back(i);
                SingValue.push_back(PosSing[TestPos]);
            }  

            Features.clear();
            for (size_t i=0;i<Faces.size();i++)
            for (size_t j=0;j<Faces[i].size();j++)
            {
                int IndexV0=Faces[i][j];
                int IndexV1=Faces[i][(j+1)%Faces[i].size()];
                Geo::Point3<ScalarType> P0=Pos[IndexV0];
                Geo::Point3<ScalarType> P1=Pos[IndexV1];
                std::pair<Geo::Point3<ScalarType>,Geo::Point3<ScalarType> > KeyE(std::min(P0,P1),std::max(P0,P1));
                if (FeaturePos.count(KeyE)==0)continue;
                std::pair<int,int> KeyF(std::min(IndexV0,IndexV1),std::max(IndexV0,IndexV1));
                Features.push_back(KeyF);
            }
            std::sort(Features.begin(), Features.end());
            Features.erase(std::unique(Features.begin(), Features.end()), Features.end());
        }

    };
}

#endif
