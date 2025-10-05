#ifndef IO_OBJ
#define IO_OBJ

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

void parseFaceLine(const std::string& line, std::vector<int>& vertexIndices) 
{
    std::stringstream ss(line);
    std::string token;
    
    while (ss >> token) {
        std::string vertexIndex;
        size_t pos = token.find('/');
        
        // Get the vertex index (before the first '/')
        if (pos != std::string::npos) {
            vertexIndex = token.substr(0, pos);
        } else {
            vertexIndex = token;
        }
        int IndexV=std::stoi(vertexIndex);
        IndexV--;
        // Convert to integer and store in the vector
        vertexIndices.push_back(IndexV);
    }
    // std::cout<<"String "<<line.c_str()<<std::endl;
    // std::cout<<"Read Face"<<std::endl;
    // for (size_t i=0;i<vertexIndices.size();i++)
    //     std::cout<<vertexIndices[i]<<","<<std::endl;
    // exit(0);
}

template <class ScalarType>
bool LoadOBJ(std::string meshFile,
             std::vector<std::vector<ScalarType > > &VertPos,
             std::vector<std::vector<ScalarType > > &VertUV,
             std::vector<std::vector<int> > &Connectivity,
             bool writeInfo=false)
{
    if (writeInfo)
        std::cout<<"Loading:"<<meshFile.c_str()<<std::endl;
    std::ifstream fin(meshFile.c_str());
    if(!fin) return false;

    //    int numV=0;
    //    int numF=0;
    std::string str;

    VertUV.clear();
    VertPos.clear();
    Connectivity.clear();

    while (!fin.eof())
    {
        std::string str;
        fin>>str;
        if (str==std::string("v"))
        {
            //numV++;
            ScalarType X,Y,Z;
            fin>>X;
            fin>>Y;
            fin>>Z;
            std::vector<ScalarType > Pos;
            Pos.push_back(X);
            Pos.push_back(Y);
            Pos.push_back(Z);
            VertPos.push_back(Pos);
        }
        if (str==std::string("vt"))
        {
            //numV++;
            ScalarType X,Y;
            fin>>X;
            fin>>Y;
            std::vector<ScalarType > Pos;
            Pos.push_back(X);
            Pos.push_back(Y);
            VertUV.push_back(Pos);
        }
        if (str==std::string("f"))
        {
            //numF++;
            std::string  line;
            std::getline(fin, line);
            std::vector<int > FIdx;
            parseFaceLine(line,FIdx); 
            Connectivity.push_back(FIdx);
        }
        //skip the line elsewhere
        if (str!=std::string("f")&&str!=std::string("v"))
        {
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }
    if (writeInfo)
    {
        std::cout<<"Read:"<<VertPos.size()<<" vertices"<<std::endl;
        std::cout<<"Read:"<<Connectivity.size()<<" faces"<<std::endl;
    }
    fin.close();
    return true;
}

template <class ScalarType>
bool LoadOBJ(std::string meshFile,
             std::vector<Geo::Point3<ScalarType> > &VertPos,
             std::vector<std::vector<int> > &Connectivity,
             bool writeInfo=false)
{
    VertPos.clear();
    Connectivity.clear();

    std::vector<std::vector<ScalarType > > VertUVVect;
    std::vector<std::vector<ScalarType > > VertPosVect;
    bool has_loaded=LoadOBJ(meshFile,VertPosVect,VertUVVect,Connectivity,writeInfo);

    if (has_loaded)
    {
        for (size_t i=0;i<VertPosVect.size();i++)
        {
            assert(VertPosVect[i].size()==3);
            VertPos.push_back(Geo::Point3<ScalarType>(VertPosVect[i][0],
                              VertPosVect[i][1],
                    VertPosVect[i][2]));
        }
    }

    return has_loaded;
}

template <class ScalarType>
bool LoadOBJ(std::string meshFile,
             std::vector<Geo::Point3<ScalarType> > &VertPos,
             std::vector<Geo::Point2<ScalarType> > &VertUV,
             std::vector<std::vector<int> > &Connectivity,
             bool writeInfo=false)
{
    VertPos.clear();
    VertUV.clear();
    Connectivity.clear();

    std::vector<std::vector<ScalarType > > VertPosVect;
    std::vector<std::vector<ScalarType > > VertUVVect;
    bool has_loaded=LoadOBJ(meshFile,VertPosVect,VertUVVect,Connectivity,writeInfo);

    if (has_loaded)
    {
        for (size_t i=0;i<VertPosVect.size();i++)
        {
            assert(VertPosVect[i].size()==3);
            VertPos.push_back(Geo::Point3<ScalarType>(VertPosVect[i][0],
                              VertPosVect[i][1],
                    VertPosVect[i][2]));
        }
        for (size_t i=0;i<VertUVVect.size();i++)
        {
            assert(VertUVVect[i].size()==2);
            VertUV.push_back(Geo::Point2<ScalarType>(VertUVVect[i][0],
                              VertUVVect[i][1]));
        }
    }

    return has_loaded;
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const std::vector<std::vector<ScalarType > > &VertPos,
              const std::vector<std::vector<int> > &Connectivity)
{
    //    std::cout<<"Loading:"<<meshFile.c_str()<<std::endl;
    std::ofstream fout(meshFile.c_str());
    if(!fout) return false;

    fout<<std::setprecision(17);

    for (size_t i=0;i<VertPos.size();i++)
    {
        assert(VertPos[i].size()==3);
        fout<<"v "<<VertPos[i][0]<<" "<<VertPos[i][1]<<" "<<VertPos[i][2]<<"\n";
    }

    for (size_t i=0;i<Connectivity.size();i++)
    {
        fout<<"f";
        for (size_t j=0;j<Connectivity[i].size();j++)
            fout<<" "<<Connectivity[i][j]+1;
        fout<<"\n";
    }

    fout.close();
    return true;
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const std::vector<std::vector<ScalarType > > &VertPos,
              const std::vector<std::vector<ScalarType > > &VertUVs,
              const std::vector<std::vector<int> > &Connectivity,
              const std::string pathTxt)
{
    //    std::cout<<"Loading:"<<meshFile.c_str()<<std::endl;
    std::ofstream fout(meshFile.c_str());
    if(!fout) return false;
    
    std::string MTLFile=RemoveExtension(meshFile);
    MTLFile+=".mtl";

    if (!pathTxt.empty())
        fout<<"mtllib "<<MTLFile.c_str()<<"\n";

    fout<<std::setprecision(17);

    for (size_t i=0;i<VertPos.size();i++)
    {
        assert(VertPos[i].size()==3);
        fout<<"v "<<VertPos[i][0]<<" "<<VertPos[i][1]<<" "<<VertPos[i][2]<<"\n";
    }

    for (size_t i=0;i<VertUVs.size();i++)
    {
        assert(VertUVs[i].size()==1);
        fout<<"vt "<<VertUVs[i][0]<<" "<<VertUVs[i][1]<<"\n";
    }

    if (!pathTxt.empty())
        fout<<"usemtl material_0\n";

    for (size_t i=0;i<Connectivity.size();i++)
    {
        std::string KeyElem("f");
        if (Connectivity[i].size()==2)
            KeyElem=std::string("l");

        fout<<KeyElem.c_str();
        for (size_t j=0;j<Connectivity[i].size();j++)
            fout<<" "<<Connectivity[i][j]+1<<"/"<<Connectivity[i][j]+1;
        fout<<"\n";
    }

    fout.close();

     if (!pathTxt.empty())
     {
        std::ofstream fout(MTLFile.c_str());
        fout<<"newmtl material_0\n";
        fout<<"Ka 0.200000 0.200000 0.200000\n";
        fout<<"Kd 0.752941 0.752941 0.752941\n";
        fout<<"Ks 1.000000 1.000000 1.000000\n";
        fout<<"Tr 0.000000\n";
        fout<<"illum 2\n";
        fout<<"Ns 0.000000\n";
        std::string TxtFile=GetFilaname(pathTxt);
        TxtFile+="\n";
        fout<<"map_Kd "<<TxtFile.c_str()<<"\n";
     }
    return true;
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const std::vector<Geo::Point3<ScalarType> > &VertPos,
              const std::vector<std::vector<int> > &Connectivity)
{
    std::vector<std::vector<ScalarType > > VertPosVect(VertPos.size());
    for (size_t i=0;i<VertPos.size();i++)
    {
        VertPosVect[i].push_back(VertPos[i].X);
        VertPosVect[i].push_back(VertPos[i].Y);
        VertPosVect[i].push_back(VertPos[i].Z);
    }
    return (WriteOBJ(meshFile,VertPosVect,Connectivity));
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const std::vector<Geo::Point3<ScalarType> > &VertPos,
              const std::vector<Geo::Point2<ScalarType> > &UVPos,
              const std::vector<std::vector<int> > &Connectivity,
              const std::string pathTxt=std::string())
{
    assert(UVPos.size()==VertPos.size());
    std::vector<std::vector<ScalarType > > VertPosVect(VertPos.size());
    std::vector<std::vector<ScalarType > > VertUVVect(UVPos.size());
    for (size_t i=0;i<VertPos.size();i++)
    {
        VertPosVect[i].push_back(VertPos[i].X);
        VertPosVect[i].push_back(VertPos[i].Y);
        VertPosVect[i].push_back(VertPos[i].Z);
    }
    for (size_t i=0;i<VertUVVect.size();i++)
    {
        VertUVVect[i].push_back(UVPos[i].X);
        VertUVVect[i].push_back(UVPos[i].Y);
    }
    return (WriteOBJ(meshFile,VertPosVect,VertUVVect,Connectivity,pathTxt));
}

template <class ScalarType>
bool WriteOBJ(std::string meshFile,
              const std::vector<Geo::Point2<ScalarType> > &VertPos,
              const std::vector<std::vector<int> > &Connectivity)
{
    std::vector<std::vector<ScalarType > > VertPosVect(VertPos.size());
    for (size_t i=0;i<VertPos.size();i++)
    {
        VertPosVect[i].push_back(VertPos[i].X);
        VertPosVect[i].push_back(VertPos[i].Y);
        VertPosVect[i].push_back(0);
    }
    return (WriteOBJ(meshFile,VertPosVect,Connectivity));
}

// template <class ScalarType>
// bool WriteOFF(std::string meshFile,
//               const std::vector<Geo::Point3<ScalarType> > &VertPos,
//               const std::vector<std::vector<int> > &Connectivity)
// {
//     std::ofstream out(meshFile);
//     if (!out.is_open())
//         return false;

//     out << "OFF\n";
//     out << VertPos.size() << " " << Connectivity.size() << " 0\n";

//     for (const auto& p : VertPos)
//         out << p.X << " " << p.Y << " " << p.Z << "\n";

//     for (const auto& face : Connectivity) {
//         out << face.size();
//         for (int idx : face)
//             out << " " << idx;
//         out << "\n";
//     }

//     out.close();
//     return true;
// }

#endif
