#ifndef INIT_CROSS_FIELD
#define INIT_CROSS_FIELD

#include "smooth_cross_field.h"

enum SmoothFieldMode{SMFieldPolyVect,SMFieldIterative};

template <class ScalarType>
struct InitCrossFieldParam
{
  SmoothFieldMode SFMode;
  int Nring;
  bool align_borders;
  ScalarType GlobalSmoothVal;

  InitCrossFieldParam()
  {
    SFMode=SMFieldPolyVect;
    Nring=3;
    align_borders = true;
    GlobalSmoothVal = 100;
  }
};

template <class ScalarType>
void InitMeshCrossField(const std::vector<Geo::Point3<ScalarType>> &TriVertPos,
const std::vector<std::vector<int>> &TriConnectivity,
const std::vector<std::pair<int, int>> &TriFeatures,
std::vector<Field::CrossF<ScalarType>> &VertCurv,
std::vector<Field::CrossF<ScalarType>> &FaceCurv,
InitCrossFieldParam<ScalarType> Param=InitCrossFieldParam<ScalarType>())
{
  //INITIALIZE WITH CURVATURE
  std::vector<Geo::Point3<ScalarType>> FaceNormals, VertNormals;
  ComputeNormals(TriVertPos, TriConnectivity, 
                 FaceNormals, VertNormals);
  ComputeCurvatureField<ScalarType>(TriVertPos, TriConnectivity,
                               FaceNormals, VertNormals, 
                               VertCurv, FaceCurv,Param.Nring);
  
  //THEN SMOOTH
  std::vector<Field::CrossF<ScalarType>> OldFCurv = FaceCurv;
  std::vector<Field::CrossF<ScalarType>> OldVCurv = VertCurv;

  FieldSmoothParam<ScalarType> SParam;
  SParam.FixedCross.clear();
  SParam.Features = TriFeatures;
  SParam.GlobalSmoothVal=Param.GlobalSmoothVal;
  SParam.align_borders=Param.align_borders;

  ComputeNormals(TriVertPos, TriConnectivity, 
                 FaceNormals, VertNormals);
  
  if (Param.SFMode==SMFieldIterative)
    SmoothCrossIterative(TriVertPos, TriConnectivity, 
                    FaceNormals, FaceCurv, SParam);

  if (Param.SFMode==SMFieldPolyVect)
    SmoothGlobalPolyvector(TriVertPos, TriConnectivity,
                           FaceNormals, FaceCurv,SParam);
  
  RedistributeAnisotropy<ScalarType>(OldVCurv, VertCurv);
  RedistributeAnisotropy<ScalarType>(OldFCurv, FaceCurv);

  InitCrossFieldQualityAsAnisotropy<ScalarType>(VertCurv, 0);
  InitCrossFieldQualityAsAnisotropy<ScalarType>(FaceCurv, 0);
}

#endif
