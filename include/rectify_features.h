#ifndef RECTIFY_FEATURES
#define RECTIFY_FEATURES

#include "point3.h"
#include "solve/least_squares_lin_constraint_solver.h"
#include "solve/linear_constraint.h"
#include "solve/linear_constraints_functions.h"
#include <Eigen/Sparse>
#include <cassert>
#include <cotangent.h>
#include <iostream>
#include <limits>
#include <mesh_functions.h>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

namespace Geo {

template <class ScalarType>
void RectifyFeatures(std::vector<Geo::Point3<ScalarType>> &VertPos,
                        const std::vector<std::vector<int>> &Faces,
                        const std::vector<std::pair<int, int>> &Features,
                        const std::vector<int> &Corners,
                        const ScalarType CircleCurvInterval=80,
                        const ScalarType PolyCurvInterval=40,
                        const bool FixCorner=true) {

  //STEP 1 INITIALIZE LAPLACIAN CONSTRAINTS
  LinearConstraintLeastSQSolver<ScalarType> LQSolver;
  std::vector<LinearConstraint<ScalarType>> LaplaceConstr;
  GetLaplacianLinearConstraint(VertPos, Faces, false, true, LaplaceConstr);
 
  //STEP 2 GET CORNER TO BE FIXED
  std::vector<int> NewCorners;
  EdgeMeshFunctions<ScalarType>::SubSampleRegularCurvature(
      VertPos, Features, Corners, 
      NewCorners, CircleCurvInterval,
      PolyCurvInterval);

  std::vector<int> AllCorners = NewCorners;
  AllCorners.insert(AllCorners.end(), Corners.begin(), Corners.end());
  std::vector<LinearConstraint<ScalarType>> LinInterpConstr;
  GetLineStraightConstraints(VertPos, Faces, Features, AllCorners, LinInterpConstr);
  
  //THEN GET THE ONE TO BE FIXED
  std::vector<LinearConstraint<ScalarType>> FixedLinConstr;
  if (FixCorner)
    GetFixedConstraint(VertPos, AllCorners, FixedLinConstr);
  else
  {
    std::vector<int> To_Fix;
    Geo::GetFixedBoxIndex(VertPos,Faces,To_Fix);
    //std::cout<<"Fixing:"<<To_Fix.size()<<std::endl;
    GetFixedConstraint(VertPos, To_Fix, FixedLinConstr);
  } 
  
  //THEN ASSEMBLE THE CONSTRAINTS
  std::vector<LinearConstraint<ScalarType>> AllLinConstr = LaplaceConstr;

  AllLinConstr.insert(AllLinConstr.end(), FixedLinConstr.begin(),
                      FixedLinConstr.end());
  AllLinConstr.insert(AllLinConstr.end(), LinInterpConstr.begin(),
                      LinInterpConstr.end());

  //INITIALIZE THE SOLVER
  LQSolver.Init(AllLinConstr);

  //FIND A SOLUTION
  std::vector<ScalarType> Solution;
  LQSolver.Solve(Solution);

  //THEN UPDATE POSITIONS
  for (size_t i = 0; i < VertPos.size(); i++) {
    int IndexX = (i * 3);
    int IndexY = IndexX + 1;
    int IndexZ = IndexX + 2;
    ScalarType ValX = Solution[IndexX];
    ScalarType ValY = Solution[IndexY];
    ScalarType ValZ = Solution[IndexZ];
    VertPos[i] = Geo::Point3<ScalarType>(ValX, ValY, ValZ);
  }
}

} // namespace Geo

#endif
