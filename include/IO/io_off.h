#ifndef IO_OFF
#define IO_OFF

#include <cstddef>
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "../point3.h"
#include "file_path_functions.h"


template <class ScalarType>
bool WriteOFF(std::string meshFile,
              const std::vector<Geo::Point3<ScalarType> > &VertPos,
              const std::vector<std::vector<int> > &Connectivity)
{
    std::ofstream out(meshFile);
    if (!out.is_open())
        return false;

    out << "OFF\n";
    out << VertPos.size() << " " << Connectivity.size() << " 0\n";

    for (const auto& p : VertPos)
        out << p.X << " " << p.Y << " " << p.Z << "\n";

    for (const auto& face : Connectivity) {
        out << face.size();
        for (int idx : face)
            out << " " << idx;
        out << "\n";
    }

    out.close();
    return true;
}

#endif
