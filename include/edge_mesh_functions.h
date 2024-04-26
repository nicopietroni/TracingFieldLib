#ifndef EDGE_MESH_FUNCTIONS
#define EDGE_MESH_FUNCTIONS

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include "point3.h"
#include "mesh_functions.h"
#include "connectivity.h"
#include <pos3.h>
#include <plane3.h>
// #include "Space_Query/triangle_grid_2D.h"

namespace Geo
{
    template <class ScalarType>
    class EdgeMeshFunctions
    {

    public:

        static void ExtractFromMeshBoundary(const std::vector<Geo::Point3<ScalarType>> &VertPos,
                                            const std::vector<std::vector<int>> &Faces,
                                            std::vector<Geo::Point3<ScalarType>> &EdgePos,
                                            std::vector<std::pair<int, int>> &Edges)
        {
            typedef std::pair<Geo::Point3<ScalarType>, Geo::Point3<ScalarType>> PosPair;
            EdgePos.clear();
            Edges.clear();
            std::vector<std::vector<int>> BorderSeq;
            getBorderSequncesVert(Faces, BorderSeq);
            std::map<Geo::Point3<ScalarType>, int> PosIDX;
            std::set<PosPair> AddedPos;
            for (size_t i = 0; i < BorderSeq.size(); i++)
            {
                int sizeB = BorderSeq[i].size();
                for (size_t j = 0; j < BorderSeq[i].size(); j++)
                {
                    int IndexV0 = BorderSeq[i][j];
                    int IndexV1 = BorderSeq[i][(j + 1) % sizeB];
                    Geo::Point3<ScalarType> Pos0 = VertPos[IndexV0];
                    Geo::Point3<ScalarType> Pos1 = VertPos[IndexV1];
                    int NewIndexV0 = -1;
                    int NewIndexV1 = -1;
                    if (PosIDX.count(Pos0) > 0)
                        NewIndexV0 = PosIDX[Pos0];
                    else
                    {
                        EdgePos.push_back(Pos0);
                        NewIndexV0 = EdgePos.size() - 1;
                    }

                    if (PosIDX.count(Pos1) > 0)
                        NewIndexV1 = PosIDX[Pos1];
                    else
                    {
                        EdgePos.push_back(Pos1);
                        NewIndexV1 = EdgePos.size() - 1;
                    }
                    assert(NewIndexV0 >= 0);
                    assert(NewIndexV1 >= 0);

                    PosPair Key(std::min(Pos0, Pos1), std::max(Pos0, Pos1));
                    if (AddedPos.count(Key)>0)continue;
                    Edges.push_back(std::pair<int,int>(NewIndexV0,NewIndexV1));
                    AddedPos.insert(Key);
                }
            }
        }
    };
};
#endif
