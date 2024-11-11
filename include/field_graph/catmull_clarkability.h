#ifndef CATMULL_CLARKABILITY
#define CATMULL_CLARKABILITY

#include <vector>
#include <iostream>

template <class ScalarType>
ScalarType OrderedCC(ScalarType L1,
                     ScalarType L2,
                     ScalarType L3)
{
    ScalarType CC=L1+L2-L3;
    //CC+=std::min(L1,std::min(L2,L3));
    //if (CC<0)return 0;
    return (CC);
}

template <class ScalarType>
ScalarType OrderedCC(ScalarType L1,
                     ScalarType L2,
                     ScalarType L3,
                     ScalarType L4,
                     ScalarType L5)
{
    ScalarType CC=(L3+L4+L5)-(L1+L2);
    //CC+=std::min(L1,std::min(L2,std::min(L3,std::min(L4,L5))));
    //if (CC<0)return 0;
    return CC;
}


template <class ScalarType>
ScalarType CatmullClarkability3(const ScalarType &L1,
                                const ScalarType &L2,
                                const ScalarType &L3,
                                const ScalarType &SideThr)
{
    ScalarType minL=std::min(L1,std::min(L2,L3));
    ScalarType CC0=OrderedCC(L1,L2,L3);
    ScalarType CC1=OrderedCC(L2,L3,L1);
    ScalarType CC2=OrderedCC(L3,L1,L2);
    ScalarType MinCC=std::min(CC0,std::min(CC1,CC2));
    MinCC+=minL*SideThr;
    //return (std::min(CC0,std::min(CC1,CC2)));
    return (MinCC);
}

template <class ScalarType>
ScalarType CatmullClarkability4(const ScalarType &L1,
                                const ScalarType &L2,
                                const ScalarType &L3,
                                const ScalarType &L4,
                                const ScalarType &SideThr)
{
    ScalarType minL=std::min(L1,std::min(L2,std::min(L3,L4)));
    ScalarType CC0=std::min(L1,L3)-std::max(L1,L3);
    ScalarType CC1=std::min(L2,L4)-std::max(L2,L4);
    ScalarType CC=std::min(CC0,CC1);
    CC+=minL*SideThr;
    //CC+=std::min(L1,std::min(L2,std::min(L3,L4)));
    return CC;

}

template <class ScalarType>
ScalarType CatmullClarkability5(const ScalarType &L1,
                                const ScalarType &L2,
                                const ScalarType &L3,
                                const ScalarType &L4,
                                const ScalarType &L5,
                                const ScalarType &SideThr)
{
    ScalarType minL=std::min(L1,std::min(L2,std::min(L3,std::min(L4,L5))));
    ScalarType CC0=OrderedCC(L1,L2,L3,L4,L5);
    ScalarType CC1=OrderedCC(L2,L3,L4,L5,L1);
    ScalarType CC2=OrderedCC(L3,L4,L5,L1,L2);
    ScalarType CC3=OrderedCC(L4,L5,L1,L2,L3);
    ScalarType CC4=OrderedCC(L5,L1,L2,L3,L4);
    ScalarType MinCC=std::min(CC0,std::min(CC1,std::min(CC2,std::min(CC3,CC4))));
    MinCC+=minL*SideThr;
    return (MinCC);
}

template <class ScalarType>
ScalarType CatmullClarkability6(const ScalarType &L1,
                                const ScalarType &L2,
                                const ScalarType &L3,
                                const ScalarType &L4,
                                const ScalarType &L5,
                                const ScalarType &L6,
                                const ScalarType &SideThr)
{
    ScalarType minL=std::min(L1,std::min(L2,std::min(L3,std::min(L4,std::min(L5,L6)))));
    ScalarType CC0=OrderedCC(L1,L3,L5);
    ScalarType CC1=OrderedCC(L2,L4,L6);
    ScalarType MinCC=std::min(CC0,CC1);
    MinCC+=minL*SideThr;
    return (std::min(CC0,CC1));
}

template <class ScalarType>
bool IsCatmullClarkable(const std::vector<ScalarType> &EdgeL,
                        const ScalarType &SideThr,
                        bool SkipValence4,
                        bool print_debug=false)
{
    if (EdgeL.size()<3)return false;
    if (EdgeL.size()>6)return false;

    if ((SkipValence4)&&(EdgeL.size()==4))return true;

    ScalarType CC=0;
    if (EdgeL.size()==3)
        CC=(CatmullClarkability3(EdgeL[0],EdgeL[1],EdgeL[2],SideThr));
    if (EdgeL.size()==4)
        CC=(CatmullClarkability4(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3],SideThr));
    if (EdgeL.size()==5)
        CC=(CatmullClarkability5(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3],EdgeL[4],SideThr));
    if (EdgeL.size()==6)
        CC=(CatmullClarkability6(EdgeL[0],EdgeL[1],EdgeL[2],EdgeL[3],EdgeL[4],EdgeL[5],SideThr));

    //CC+=SideThr;

    //  if (CC<0)
    // {
    //     std::cout<<"*** CC Computation ***"<<CC<<std::endl;
    //     std::cout<<"Side "<<EdgeL.size()<<std::endl;
    //     for (size_t i=0;i<EdgeL.size();i++)
    //         std::cout<<EdgeL[i]<<",";
    //     std::cout<<std::endl;

    //     std::cout<<"SideThr: "<<SideThr<<std::endl;
    //     std::cout<<"CC: "<<CC<<std::endl;
    // }
    //CC=std::max(CC,(ScalarType)0);
    //    if ((EdgeL.size()==4)&&(CC>Ratio4))
    //        return true;
    //    if ((EdgeL.size()==4)&&(CC<=Ratio4))
    //        return false;

    return (CC>0);
}

//template <class ScalarType>
//size_t AddedSingularities(//const int NumF,
//                          const std::vector<ScalarType> &EdgeL,
//                          const ScalarType &SideThr,
//                          bool SkipValence4,
//                          bool print_debug=false)
//{
//    assert(SideThr>0);
//    //bool CC=IsCatmullClarkable(NumF,EdgeL,SideThr,SkipValence4,print_debug);
//    bool CC=IsCatmullClarkable(EdgeL,SideThr,SkipValence4,print_debug);
//    //std::cout<<"test CC "<<CC<<std::endl;
//    if (CC)
//        return 0;
//    else
//        return 1;
//}

#endif
