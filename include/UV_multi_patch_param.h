#ifndef UV_MULTI_PATCH_PARAM
#define UV_MULTI_PATCH_PARAM

#include <point2.h>
#include <point3.h>
#include <UV_param.h>
#include <mesh_functions.h>
#include <outline_packer.h>
#include <IO/io_obj.h>

namespace Geo
{

    template <class ScalarType>
    class UVMultiPatchParam
    {

    public:
        // enum UVType
        // {
        //     UVArap,
        //     UVLSConformal,
        //     UVTutte
        // };

        struct MultiPatchParameters
        {
            ScalarType border_size;
            bool rearrange;
            int Max_iteNum; // for Arap
            typename Geo::UVUtils<ScalarType>::UVType ParamType;
            std::vector<int> CornerIdx; // for Tutte
            bool solve_degenerate;

            MultiPatchParameters()
            {
                border_size = 5;
                rearrange = true;
                Max_iteNum = 100;
                ParamType = Geo::UVUtils<ScalarType>::UVLSConformal;
                solve_degenerate = true;
            }
        };

        static bool Parametrize(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                const std::vector<std::vector<int>> &Faces,
                                std::vector<Geo::Point2<ScalarType>> &VertUV,
                                MultiPatchParameters MPatchParam = MultiPatchParameters())
        {
            std::vector<std::vector<Geo::Point3<ScalarType>>> SubMeshesV;
            std::vector<std::vector<std::vector<int>>> SubMeshesF;
            std::vector<std::vector<int>> OriginalV, OriginalF;
            SplitByFFConnectedComponenent<ScalarType>(VertPos, Faces, SubMeshesV,
                                                      SubMeshesF, OriginalV, OriginalF);

            assert(SubMeshesF.size() == SubMeshesV.size());
            assert(SubMeshesF.size() == OriginalV.size());

            VertUV.clear();
            VertUV.resize(VertPos.size(), Geo::Point2<ScalarType>(0, 0));

            std::vector<Field::CrossF<ScalarType>> FaceCurv;
            for (size_t i = 0; i < SubMeshesV.size(); i++)
            {
                std::vector<Geo::Point2<ScalarType>> CompVertUV;
                bool solved = false;
                // if (MPatchParam.solve_degenerate)
                //     SmoothDegerateTriangles<ScalarType>(SubMeshesV[i], SubMeshesF[i], 1);

                switch (MPatchParam.ParamType)
                {
                case Geo::UVUtils<ScalarType>::UVLSConformal:

                    solved = Geo::UVParam<ScalarType>::LeastSquaresConformal(SubMeshesV[i],
                                                                             SubMeshesF[i],
                                                                             CompVertUV);
                    {
                        Geo::Box2<ScalarType> BBUV;
                        BBUV.Init(CompVertUV);
                        if (std::isnan(BBUV.DimX()))
                        {
                            std::cout<<"SOLVING"<<std::endl;
                            //SmoothLaplacian<ScalarType>(SubMeshesV[i], SubMeshesF[i], 1);
                            //SmoothDegerateTriangles<ScalarType>(SubMeshesV[i], SubMeshesF[i], 1);
                            Geo::UVParam<ScalarType>::LeastSquaresConformal(SubMeshesV[i],
                                                                            SubMeshesF[i],
                                                                            CompVertUV);
                            Geo::Box2<ScalarType> BBUV;
                            BBUV.Init(CompVertUV);
                            if (std::isnan(BBUV.DimX()))
                            {
                                WriteOBJ("testNaN.obj", SubMeshesV[i], SubMeshesF[i]);
                                std::cout << "DE boia" << std::endl;
                                //exit(0);
                            }
                        }
                    }
                    break;
                case Geo::UVUtils<ScalarType>::UVArap:
                    solved = Geo::UVParam<ScalarType>::Arap(SubMeshesV[i], SubMeshesF[i], CompVertUV,
                                                            true, MPatchParam.Max_iteNum);
                    // if (!solved)
                    //     return false;
                    break;

                case Geo::UVUtils<ScalarType>::UVPoisson:

                    solved = Geo::UVParam<ScalarType>::Poisson(SubMeshesV[i], SubMeshesF[i], FaceCurv,
                                                               false, 1, CompVertUV);
                    break;
                default:

                    assert(MPatchParam.ParamType == Geo::UVUtils<ScalarType>::UVTutte);
                    std::vector<std::vector<int>> SubCorner;
                    std::set<int> CornerIdxSet(MPatchParam.CornerIdx.begin(),
                                               MPatchParam.CornerIdx.end());

                    SubCorner.resize(OriginalV.size());

                    for (size_t i = 0; i < OriginalV.size(); i++)
                        for (size_t j = 0; j < OriginalV[i].size(); j++)
                        {
                            int OrIndex = OriginalV[i][j];
                            if (CornerIdxSet.count(OrIndex) == 0)
                                continue;
                            SubCorner[i].push_back(j);
                        }
                    std::vector<int> CornerIDX = SubCorner[i];
                    solved = Geo::UVParam<ScalarType>::TutteEmbedding(SubMeshesV[i], SubMeshesF[i],
                                                                      CornerIDX, CompVertUV, false);
                }
                if (!solved)
                    return false;
                // then reassign
                assert(OriginalV[i].size() == SubMeshesV[i].size());
                assert(OriginalV[i].size() == CompVertUV.size());
                for (size_t j = 0; j < OriginalV[i].size(); j++)
                {
                    int GlobalIdx = OriginalV[i][j];
                    assert(GlobalIdx >= 0);
                    assert(GlobalIdx < VertUV.size());
                    VertUV[GlobalIdx] = CompVertUV[j];
                }
            }
            Geo::Box2<ScalarType> BB_UV;
            BB_UV.Init(VertUV);
            if (MPatchParam.rearrange)
            {
                Geo::OutlinePacker<ScalarType> Packer;
                Packer.Pack(VertUV, Faces, MPatchParam.border_size);
            }
            return true;
        }

        // static bool LeastSquaresConformal(const std::vector<Geo::Point3<ScalarType> > &VertPos,
        //                                   const std::vector<std::vector<int> > &Faces,
        //                                   std::vector<Geo::Point2<ScalarType> > &VertUV,
        //                                   bool Rearrange=true)
        // {
        //     std::vector<std::vector<Geo::Point3<ScalarType> > > SubMeshesV;
        //     std::vector<std::vector<std::vector<int> > > SubMeshesF;
        //     std::vector<std::vector<int> > OriginalV,OriginalF;
        //     SplitByFFConnectedComponenent<ScalarType>(VertPos,Faces,SubMeshesV,
        //                                               SubMeshesF,OriginalV,OriginalF);

        //     assert(SubMeshesF.size()==SubMeshesV.size());
        //     assert(SubMeshesF.size()==OriginalV.size());

        //     VertUV.clear();
        //     VertUV.resize(VertPos.size(),Geo::Point2<ScalarType>(0,0));
        //     for (size_t i=0;i<SubMeshesV.size();i++)
        //     {
        //         std::vector<Geo::Point2<ScalarType> > CompVertUV;
        //         //std::cout<<"Parametrizing Partition "<<i<<std::endl;

        //         //WriteOBJ<ScalarType>(std::string("test.obj"),SubMeshesV[i],SubMeshesF[i]);
        //         bool solved=Geo::UVParam<ScalarType>::LeastSquaresConformal(SubMeshesV[i],SubMeshesF[i],CompVertUV);
        //         if (!solved)return false;

        //         //then reassign
        //         assert(OriginalV[i].size()==SubMeshesV[i].size());
        //         assert(OriginalV[i].size()==CompVertUV.size());
        //         for (size_t j=0;j<OriginalV[i].size();j++)
        //         {
        //             int GlobalIdx=OriginalV[i][j];
        //             assert(GlobalIdx>=0);
        //             assert(GlobalIdx<VertUV.size());
        //             VertUV[GlobalIdx]=CompVertUV[j];
        //         }
        //     }
        //     Geo::Box2<ScalarType> BB_UV;
        //     BB_UV.Init(VertUV);
        //     if (Rearrange)
        //     {
        //         Geo::OutlinePacker<ScalarType> Packer;
        //         Packer.Pack(VertUV,Faces,BB_UV.Diag()*0.01);
        //     }

        //     return true;
        // }

        //     static bool Arap(const std::vector<Geo::Point3<ScalarType> > &VertPos,
        //                      const std::vector<std::vector<int> > &Faces,
        //                      std::vector<Geo::Point2<ScalarType> > &VertUV,
        //                      bool Rearrange=true,bool InitializeConformal=true,
        //                      int MaxIte=100)
        //     {
        //         std::vector<std::vector<Geo::Point3<ScalarType> > > SubMeshesV;
        //         std::vector<std::vector<std::vector<int> > > SubMeshesF;
        //         std::vector<std::vector<int> > OriginalV,OriginalF;
        //         SplitByFFConnectedComponenent<ScalarType>(VertPos,Faces,SubMeshesV,
        //                                                   SubMeshesF,OriginalV,OriginalF);

        //         assert(SubMeshesF.size()==SubMeshesV.size());
        //         assert(SubMeshesF.size()==OriginalV.size());

        //         VertUV.clear();
        //         VertUV.resize(VertPos.size(),Geo::Point2<ScalarType>(0,0));
        //         for (size_t i=0;i<SubMeshesV.size();i++)
        //         {
        //             std::vector<Geo::Point2<ScalarType> > CompVertUV;
        //             bool solved=Geo::UVParam<ScalarType>::Arap(SubMeshesV[i],SubMeshesF[i],CompVertUV,InitializeConformal,MaxIte);
        //             if (!solved)return false;
        //             Geo::Box2<ScalarType> BB;
        //             BB.Init(CompVertUV);
        //             if (BB.Area()==0)
        //                 solved=Geo::UVParam<ScalarType>::LeastSquaresConformal(SubMeshesV[i],SubMeshesF[i],CompVertUV);

        //             if (!solved)return false;

        //             UVUtils<ScalarType>::ScaleToMatchArea(SubMeshesV[i],SubMeshesF[i],CompVertUV);

        //                 //WriteOBJ<ScalarType>(std::string("testArap.obj"),SubMeshesV[i],SubMeshesF[i]);
        //             //then reassign
        //             assert(OriginalV[i].size()==SubMeshesV[i].size());
        //             assert(OriginalV[i].size()==CompVertUV.size());
        //             for (size_t j=0;j<OriginalV[i].size();j++)
        //             {
        //                 int GlobalIdx=OriginalV[i][j];
        //                 assert(GlobalIdx>=0);
        //                 assert(GlobalIdx<VertUV.size());
        //                 VertUV[GlobalIdx]=CompVertUV[j];
        //             }
        //         }
        //         Geo::Box2<ScalarType> BB_UV;
        //         BB_UV.Init(VertUV);
        //         if (Rearrange)
        //         {
        //             Geo::OutlinePacker<ScalarType> Packer;
        //             Packer.Pack(VertUV,Faces,BB_UV.Diag()*0.01);
        //         }
        // //        Geo::Box2<ScalarType> BB_UV;
        // //        BB_UV.Init(VertUV);
        // //        if (Rearrange)
        // //            Geo::OutlinePacker<ScalarType>::Pack(VertUV,Faces,BB_UV.Diag()*0.01);

        //         return true;
        //     }

        //     static bool Tutte(const std::vector<Geo::Point3<ScalarType> > &VertPos,
        //                       const std::vector<std::vector<int> > &Faces,
        //                       const std::vector<int> CornerIdx,
        //                       std::vector<Geo::Point2<ScalarType> > &VertUV,
        //                       bool Rearrange=true)
        //     {
        //         std::vector<std::vector<Geo::Point3<ScalarType> > > SubMeshesV;
        //         std::vector<std::vector<std::vector<int> > > SubMeshesF;
        //         std::vector<std::vector<int> > OriginalV,OriginalF;
        //         std::vector<std::vector<int> > SubCorner;
        //         SplitByFFConnectedComponenent<ScalarType>(VertPos,Faces,SubMeshesV,
        //                                                   SubMeshesF,OriginalV,OriginalF);

        //         std::set<int> CornerIdxSet(CornerIdx.begin(),CornerIdx.end());
        //         SubCorner.resize(OriginalV.size());

        //         for (size_t i=0;i<OriginalV.size();i++)
        //             for (size_t j=0;j<OriginalV[i].size();j++)
        //             {
        //                 int OrIndex=OriginalV[i][j];
        //                 if (CornerIdxSet.count(OrIndex)==0)continue;
        //                 SubCorner[i].push_back(j);
        //             }

        //         assert(SubMeshesF.size()==SubMeshesV.size());
        //         assert(SubMeshesF.size()==OriginalV.size());

        //         VertUV.clear();
        //         VertUV.resize(VertPos.size(),Geo::Point2<ScalarType>(0,0));
        //         for (size_t i=0;i<SubMeshesV.size();i++)
        //         {
        //             std::vector<Geo::Point2<ScalarType> > CompVertUV;
        //             std::vector<int> CornerIDX=SubCorner[i];
        //             bool solved=Geo::UVParam<ScalarType>::TutteEmbedding(SubMeshesV[i],SubMeshesF[i],CornerIDX,CompVertUV);
        //             if (!solved)return false;

        //             //then reassign
        //             assert(OriginalV[i].size()==SubMeshesV[i].size());
        //             assert(OriginalV[i].size()==CompVertUV.size());
        //             for (size_t j=0;j<OriginalV[i].size();j++)
        //             {
        //                 int GlobalIdx=OriginalV[i][j];
        //                 assert(GlobalIdx>=0);
        //                 assert(GlobalIdx<VertUV.size());
        //                 VertUV[GlobalIdx]=CompVertUV[j];
        //             }
        //         }
        //         Geo::Box2<ScalarType> BB_UV;
        //         BB_UV.Init(VertUV);
        //         if (Rearrange)
        //         {
        //             Geo::OutlinePacker<ScalarType> Packer;
        //             Packer.Pack(VertUV,Faces,BB_UV.Diag()*0.01);
        //         }
        //         return true;
        //     }
    };

}

#endif
