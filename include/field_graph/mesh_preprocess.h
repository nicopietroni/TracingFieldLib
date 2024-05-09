#ifndef TRACING_MESH_PREPROCESS
#define TRACING_MESH_PREPROCESS

#include <connectivity.h>
#include <set>
#include <mesh_refine.h>

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

        static bool RefineCloseSingFacesStep(std::vector<Geo::Point3<ScalarType>> &Pos,
                                             std::vector<std::vector<int>> &Faces,
                                             std::vector<Field::CrossF<ScalarType>> &FaceCurv,
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

            std::vector<ScalarType> Alpha(To_splitE.size(), 0.5);
            std::vector<int> OldF;
            RefineMesh<ScalarType>::SplitTriangleFaces(Pos, Faces, To_splitE, Alpha, &OldF);
            assert(OldF.size() == Faces.size());

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
                                         const std::vector<int> &SingIndex,
                                         const int extend_step)
        {
            while (RefineCloseSingFacesStep(Pos, Faces, FaceCurv, SingIndex, extend_step))
            {
            }
        }

        static void DesambiguateInitialDarts(std::vector<Geo::Point3<ScalarType>> &Pos,
                                             std::vector<std::vector<int>> &Faces)
        {
            ScalarType delta= AvgEdgeLen(Pos,Faces)*0.0001;
            //ScalarType delta=std::numeric_limits<float>::epsilon()*(ScalarType)2;
            DesambiguateSamePosVertices<ScalarType>(Pos, Faces,delta);
        }
    };
}

#endif
