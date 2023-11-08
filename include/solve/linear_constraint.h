#ifndef LINEAR_CONSTRAINT
#define LINEAR_CONSTRAINT

#include <cstddef>
#include <vector>
#include <point3.h>

namespace Geo
{

template <class ScalarType>
struct VarTerm
{
    size_t IndexVar;
    size_t IndexCoord;
    ScalarType W;

    VarTerm(size_t _IndexVar,
            size_t _IndexCoord,
            ScalarType _W)
    {
        IndexVar=_IndexVar;
        IndexCoord=_IndexCoord;
        W=_W;
    }

    VarTerm()
    {
        IndexVar=0;
        IndexCoord=0;
        W=0;
    }
};

template <class ScalarType>
class LinearConstraint
{

public:

    ScalarType RHSide;
    std::vector<VarTerm<ScalarType> > LHSide;
    ScalarType Weight;

    LinearConstraint()
    {
        RHSide=0;
        Weight=1;
    }

};

//template <class ScalarType>
//void FixedLinearConstraints(const std::vector<Geo::Point3<ScalarType> > &Pos,
//                            const std::vector<int> &FixedV,
//                            const ScalarType Weight,
//                            std::vector<LinearConstraint<ScalarType> > &FixConstr)
//{
//    FixConstr.clear();
//    for (size_t i=0;i<Pos.size();i++)
//        for (size_t j=0;j<3;j++)
//        {
//            Geo::LinearConstraint<ScalarType> Constr;
//            Constr.RHSide=Pos[i].cV(j);
//            Constr.Weight=Weight;
//            Geo::VarTerm<ScalarType> Var(i,j,1);
//            Constr.LHSide.push_back(Var);
//            FixConstr.push_back(Constr);
//        }
//}

}

#endif
