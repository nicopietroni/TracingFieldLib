#ifndef IO_POLY
#define IO_POLY

#include <cassert>
#include <cstddef>
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "file_path_functions.h"
#include "../point3.h"
#include "../tetra3.h"

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

    bool foundTet=false;
    do {
    fin>>str;
    if (str==std::string("Tetrahedra"))foundTet=true;
    }while (str!=std::string("Tetrahedra"));
    if (!foundTet)return false;
    
    fin>>n;
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

        Tetra::Tetra3<ScalarType> T(VertPos[Connectivity[i][0]],
                                    VertPos[Connectivity[i][1]],
                                    VertPos[Connectivity[i][2]],
                                    VertPos[Connectivity[i][3]]);

        if(T.Volume()<0)
        {
            std::swap(Connectivity[i][0],Connectivity[i][1]);
            Tetra::Tetra3<ScalarType> T1(VertPos[Connectivity[i][0]],
                                    VertPos[Connectivity[i][1]],
                                    VertPos[Connectivity[i][2]],
                                    VertPos[Connectivity[i][3]]);
            assert(T1.Volume()>=0);
        }
       
        for (size_t j=0;j<4;j++)
        {
            assert(Connectivity[i][j]>=0);
            assert(Connectivity[i][j]<VertPos.size());
        }
        //check volume orientation

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
