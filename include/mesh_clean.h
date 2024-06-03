#ifndef MESH_CLEAN
#define MESH_CLEAN

#include "point3.h"
#include "mesh_functions.h"

namespace Geo
{
    void GetDegenerateFaces(const std::vector<std::vector<int>> &Faces,
                            std::vector<int> &DegenerateIdx)
    {
        DegenerateIdx.clear();
        for (size_t i = 0; i < Faces.size(); i++)
        {
            bool found = false;
            for (size_t j = 0; j < Faces[i].size() - 1; j++)
            {
                for (size_t k = (j + 1); k < Faces[i].size(); k++)
                {
                    if (Faces[i][j] != Faces[i][k])
                        continue;
                    DegenerateIdx.push_back(i);
                    break;
                }
                if (found)
                    break;
            }
        }
    }

    template <class ScalarType>
    void RemoveDegenerateFaces(std::vector<Geo::Point3<ScalarType>> &Pos,
                               std::vector<std::vector<int>> &Faces)
    {
        std::vector<int> DegenerateIdx;
        GetDegenerateFaces(Faces,DegenerateIdx);
        if (DegenerateIdx.size()==0)return;
        std::set<int> DegenerateIdxSet(DegenerateIdx.begin(),DegenerateIdx.end());
        std::vector<std::vector<int>> New_Faces;
        for (size_t i=0;i<Faces.size();i++)
        {
            if (DegenerateIdxSet.count(i)>0)continue;
            New_Faces.push_back(Faces[i]);
        }
        Faces=New_Faces;
        RemoveUnfererencedVert(Pos,Faces);
    }
}

#endif
