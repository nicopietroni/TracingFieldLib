#ifndef IO_POLY
#define IO_POLY

#include <cstddef>
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <point3.h>
#include <IO/file_path_functions.h>

template <class ScalarType>
bool LoadTetraMesh(std::string meshFile,
                   std::vector<Geo::Point3<ScalarType> > &VertPos,
                   std::vector<std::vector<int> > &Connectivity,
                   bool writeInfo=true)
{
    std::ifstream fin(meshFile.c_str());
    if(!fin) return false;

    int n;
    std::string str;

    fin>>str;
    fin>>n;
    assert(str==std::string("MeshVersionFormatted"));
    assert(n==1);
    fin>>str;
    fin>>n;
    assert(str==std::string("Dimension"));
    assert(n==3);
    fin>>str;
    fin>>n;
    assert(str==std::string("Vertices"));

    VertPos.clear();
    Connectivity.clear();

    //then allocate the vertices
    VertPos.resize(n);
    for (int i=0;i<n;i++)
    {
        double w;
        fin>>VertPos[i].X;
        fin>>VertPos[i].Y;
        fin>>VertPos[i].Z;
        fin>>w;
    }
    fin>>str;
    fin>>n;
    assert(str==std::string("Tetrahedra"));
    Connectivity.resize(n,std::vector<int>(4,-1));
    for (int i=0;i<n;i++)
    {
        int IndexT;
        int IndexV;
        fin>>IndexV;
        Connectivity[i][0]=IndexV-1;
        fin>>IndexV;
        Connectivity[i][1]=IndexV-1;
        fin>>IndexV;
        Connectivity[i][2]=IndexV-1;
        fin>>IndexV;
        Connectivity[i][3]=IndexV-1;
        fin>>IndexT;
        for (size_t j=0;j<4;j++)
        {
            assert(Connectivity[i][j]>=0);
            assert(Connectivity[i][j]<VertPos.size());
        }
    }
    fin.close();

    if (writeInfo)
    {
        std::cout<<"Read:"<<Connectivity.size()<<" tets"<<std::endl;
        std::cout<<"Read:"<<VertPos.size()<<" vertices"<<std::endl;
    }
    return true;
}

#endif
