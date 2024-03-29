#ifndef MESH_CREATE
#define MESH_CREATE

#include <vector>
#include "point3.h"
#include "mesh_functions.h"

// #include "Space_Query/triangle_grid_2D.h"
namespace Geo
{
    template <class ScalarType>
    class MeshCreate
    {
    public:
        static bool UnitSphericalDomeAngle(const ScalarType &angle_interval,
                                           std::vector<Geo::Point3<ScalarType>> &VertPos,
                                           std::vector<std::vector<int>> &Faces,
                                           int Subd0 = 48, int Subd1 = 48)
        {
            VertPos.clear();
            Faces.clear();
            if (angle_interval <= 0)
                return false;
            if (angle_interval >= M_PI / 2)
                return false;

            ScalarType Phi_Interval = angle_interval; // 2*(M_PI/2-angle_interval);
            ScalarType Phi_Step = Phi_Interval / Subd0;
            for (int i = 0; i <= Subd0; ++i)
            {
                // double phi = (M_PI/2 * i) / (Subd0 - 1); // Range of phi: [0, PI]
                // double phi = angle_interval+Phi_Step*i;
                double phi = Phi_Step * i;
                for (int j = 0; j < Subd1; ++j)
                {
                    double theta = 2 * M_PI * j / Subd1; // Range of theta: [0, 2PI]

                    ScalarType XVal = sin(phi) * cos(theta);
                    ScalarType YVal = sin(phi) * sin(theta);
                    ScalarType Zval = cos(phi);
                    VertPos.push_back(Geo::Point3<ScalarType>(XVal, YVal, Zval));
                }
            }
            // then find the tessellation
            for (int i = 1; i < Subd0; ++i)
                for (int j = 0; j < Subd1; ++j)
                {
                    int IndexV0 = i * Subd1 + j;
                    int IndexV1 = i * Subd1 + (j + 1) % Subd1;
                    int IndexV2 = (i + 1) * Subd1 + j;
                    int IndexV3 = (i + 1) * Subd1 + (j + 1) % Subd1;
                    Faces.resize(Faces.size() + 1);
                    Faces.back().push_back(IndexV0);
                    Faces.back().push_back(IndexV1);
                    Faces.back().push_back(IndexV2);
                    Faces.resize(Faces.size() + 1);
                    Faces.back().push_back(IndexV1);
                    Faces.back().push_back(IndexV3);
                    Faces.back().push_back(IndexV2);
                }
            // add the last level
            int LastV = 0;
            for (int j = 0; j < Subd1; ++j)
            {
                int IndexV0 = LastV;
                int IndexV1 = Subd1 + j;
                int IndexV2 = Subd1 + (j + 1) % Subd1;
                Faces.resize(Faces.size() + 1);
                Faces.back().push_back(IndexV0);
                Faces.back().push_back(IndexV2);
                Faces.back().push_back(IndexV1);
            }
            RemoveUnfererencedVert<ScalarType>(VertPos,Faces);
            return true;
        }

        static bool UnitSphericalDomeDist(const ScalarType &distance,
                                          std::vector<Geo::Point3<ScalarType>> &VertPos,
                                          std::vector<std::vector<int>> &Faces,
                                          int Subd0 = 48, int Subd1 = 48)
        {
            if (distance <= 0)
                return false;
            if (distance >= 1)
                return false;

            ScalarType angle_interval = M_PI / 2 - asin(1 - distance);
            UnitSphericalDomeAngle(angle_interval, VertPos, Faces, Subd0, Subd1);
            return true;
        }

        static bool UnitSphericalDomeDist(const ScalarType &radius_sphere,
                                          const ScalarType &distance,
                                          std::vector<Geo::Point3<ScalarType>> &VertPos,
                                          std::vector<std::vector<int>> &Faces,
                                          int Subd0 = 48, int Subd1 = 48)
        {
            if (distance <= 0)
                return false;
            if (distance >= radius_sphere)
                return false;

            ScalarType scaled_dist = distance / radius_sphere;
            bool IsOK = UnitSphericalDomeDist(scaled_dist, VertPos, Faces, Subd0, Subd1);
            for (size_t i = 0; i < VertPos.size(); i++)
                VertPos[i] = VertPos[i] * radius_sphere;

            return IsOK;
        }

        static bool SphericalDomeDist(const ScalarType &base_radius, const ScalarType &dome_height,
                                      std::vector<Geo::Point3<ScalarType>> &VertPos,
                                      std::vector<std::vector<int>> &Faces, int subSteps = 24)
        {
            // ScalarType RadiusSphere=sqrt(2 * (base_radius * base_radius + base_radius * dome_height) + dome_height * dome_height);
            // ScalarType RadiusSphere=dome_height / 2 + sqrt(pow(dome_height, 2) / 4 - pow(base_radius, 2));
            ScalarType RadiusSphere = (pow(base_radius, 2) + pow(dome_height, 2)) / (2 * dome_height);
            std::cout << "Radius:" << RadiusSphere << std::endl;
            std::cout << "Height:" << dome_height << std::endl;
            return UnitSphericalDomeDist(RadiusSphere, dome_height, VertPos, Faces, subSteps, subSteps*2);
        }
    };
}

#endif
