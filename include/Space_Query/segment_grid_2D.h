#ifndef SEGMENT_GRID_2D
#define SEGMENT_GRID_2D

#include <vector>
#include "grid_basics_2D.h"
#include <point2.h>
#include <point3.h>
#include <box2.h>
#include <segment2.h>

namespace  Geo{

/*! \brief The cell of the Grid that keep the triangles */
template <class SegmentType>
class SegmentCell2D
{
public:
    /*! \brief The reference to the primitives*/
    std::vector<int> entries;

    /*! \brief Base constructor */
    SegmentCell2D(){}

    /*! \brief Insert a primitive in the Cell
     *  \param tri The primitive to be inserted
     */
    void Insert(int tri)
    {
        entries.push_back(tri);
    }

};

/*! \brief This class implement a Grid-like search structure for triangular faces */
template <class ScalarType>
class SegmentSpatialIndex2D
{

    typedef SegmentCell2D<ScalarType> CellType;
    typedef Geo::Point2<int> GridIndex;
    typedef Geo::Point2<ScalarType> CoordType;

    /*! \brief The actual grid of cells*/
    std::vector<std::vector<CellType> > Grid;

public:
     /*! \brief The bounding box of the whole space*/
    //AABBox<ScalarType> bbox;
    Box2<ScalarType> bbox;

    /*! \brief Number of cells forming the grid*/
    GridIndex siz;
    /*! \brief Dimensions of a single cell*/
    CoordType voxel;
    /*! \brief Global temporary Mark */
    int TMark;
    /*! \brief Per segment temporary Mark */
    std::vector<int> SegMark;

private:

    /*! \brief Return the index of the cell for a given 3D point
     *  \param p The input 3D point
     *  \param pi The grid coordinate
     */
    inline void PToIP(const CoordType & p, GridIndex &pi ) const
    {
        CoordType t = p - bbox.Min;
        pi.X = int( t.X / voxel.X );
        pi.Y = int( t.Y / voxel.Y );
    }

    /*! \brief Return the primitives within a cell
     *  \param Index The grid coordinates
     *  \param ret The output primitives
     */
    inline void Get(const GridIndex &Index,
                    std::vector<int> &ret)const
    {
        ret=Grid[Index.X][Index.Y].entries;
    }

    /*! \brief Convert a Bounding box to Grid Coordinates
     *  \param b The input bounding box
     *  \param ib The output Grid bounding box
     */
    inline void BoxToIBox( const Box2<ScalarType> & b,
                           Box2<int> & ib ) const
    {
        PToIP(b.Min, ib.Min);
        PToIP(b.Max, ib.Max);
    }

    /*! \brief Insert a segment in a Given Cell
     *  \param tris The index of the triangle
     *  \param Index The Grid Coordinates
     */
    void Insert(int seg,
                GridIndex Index)
    {
        Grid[Index.X][Index.Y].Insert(seg);
    }

    /*! \brief Insert a triangle in its cells
     *  \param vert The vertex vector of the input mesh
     *  \param face The face vector of the input mesh
     *  \param IndexF The index of the face
     */
    void Insert(const std::vector<Geo::Segment2<ScalarType> > &segs,
                int IndexS)
    {
        //initialize temporary mark
        SegMark[IndexS]=TMark;

        Box2<ScalarType> segBB=SegmentBBox(segs[IndexS]);
        Box2<int> intBox;
        BoxToIBox(segBB,intBox);
        for (int X=intBox.Min.X;X<=intBox.Max.X;X++)
            for (int Y=intBox.Min.Y;Y<=intBox.Max.Y;Y++)
                    Insert(IndexS,GridIndex(X,Y));
    }

public:


    void SegmentsInBox(const std::vector<Geo::Segment2<ScalarType> > &segs,
                       const Box2<ScalarType> &testBB,
                       std::vector<int> &SegIn)const
    {
        SegIn.clear();

        Box2<ScalarType> bb_filtered=testBB;
        bb_filtered.Intersect(bbox);
        if (bb_filtered.IsEmpty())return;

        Box2<int> intBox;
        BoxToIBox(bb_filtered,intBox);

        for (int xi=intBox.Min.X;xi<=intBox.Max.X;xi++)
            for (int yi=intBox.Min.Y;yi<=intBox.Max.Y;yi++)
            {
                if(xi>=Grid.size())continue;
                if(yi>=Grid[xi].size())continue;
                SegIn.insert(SegIn.end(),
                             Grid[xi][yi].entries.begin(),
                             Grid[xi][yi].entries.end());
            }
        std::sort(SegIn.begin(), SegIn.end());
        std::vector<int>::iterator last = std::unique(SegIn.begin(), SegIn.end());
        SegIn.erase(last, SegIn.end());
    }

    void InteresectSeg(const std::vector<Geo::Segment2<ScalarType> > &segs,
                       const Geo::Segment2<ScalarType> &test_seg,
                       std::vector<int> &IntSegIdx,
                       std::vector<Geo::Point2<ScalarType> > &IntSegPos)const
    {
        IntSegIdx.clear();
        IntSegPos.clear();

        Box2<ScalarType> bb_seg=SegmentBBox(test_seg);
        std::vector<int> SegIn;
        SegmentsInBox(segs,bb_seg,SegIn);

        for (size_t i=0;i<SegIn.size();i++)
        {
            int IndexTest=SegIn[i];
            assert(IndexTest>=0);
            assert(IndexTest<segs.size());
            Geo::Segment2<ScalarType> S0=segs[IndexTest];
            Geo::Segment2<ScalarType> S1=test_seg;
            ScalarType lambda0;
            if (Geo::Intersect<ScalarType>(S0,S1,lambda0))
            {
                Geo::Point2<ScalarType> IntP=S0.P(0)*(lambda0)+S0.P(1)*(1-lambda0);
                IntSegPos.push_back(IntP);
                IntSegIdx.push_back(IndexTest);
            }
        }
    }

    bool FirstInteresect(const std::vector<Geo::Segment2<ScalarType> > &segs,
                        const Geo::Point2<ScalarType> &Orig,
                        const Geo::Point2<ScalarType> &Dir,
                        int &FirstIntSegIdx,Geo::Point2<ScalarType> &FirstIntSegPos)const
    {
        Geo::Segment2<ScalarType> test_seg(Orig,Orig+Dir);
        std::vector<int> IntSegIdx;
        std::vector<Geo::Point2<ScalarType> > IntSegPos;
        InteresectSeg(segs,test_seg,IntSegIdx,IntSegPos);
        if (IntSegIdx.size()==0)return false;
        if (IntSegIdx.size()==1)
        {
            FirstIntSegIdx=IntSegIdx[0];
            FirstIntSegPos=IntSegPos[0];
            return true;
        }

        ScalarType minD=std::numeric_limits<ScalarType>::max();
        for (size_t i=0;i<IntSegIdx.size();i++)
        {
            ScalarType durrD=(IntSegPos[i]-test_seg.P(0)).Norm();
            if (durrD>minD)continue;
            minD=durrD;
            FirstIntSegIdx=IntSegIdx[i];
            FirstIntSegPos=IntSegPos[i];
        }
        return true;
    }

    bool FirstInteresect(const std::vector<Geo::Segment2<ScalarType> > &segs,
                        const Geo::Point2<ScalarType> &Orig,
                        const Geo::Point2<ScalarType> &Dir,
                        int &FirstIntSegIdx,ScalarType &lambda)const
    {
        Geo::Point2<ScalarType> FirstIntSegPos;
        bool has_int=FirstInteresect(segs,Orig,Dir,FirstIntSegIdx,FirstIntSegPos);
        if (has_int)
        {
            assert(FirstIntSegIdx>=0);
            assert(FirstIntSegIdx<segs.size());
            Geo::Segment2<ScalarType> curr_s=segs[FirstIntSegIdx];
            ScalarType L=Lenght(curr_s);
            ScalarType D0=(FirstIntSegPos-curr_s.P(0)).Norm();
            assert(D0<=L);
            lambda=1-D0/L;
            return true;
        }
        return false;
    }

    /*! \brief Initialize the search structure
     *  \param vert The vertex vector of the input mesh
     *  \param face The face vector of the input mesh
     */
    void Init(const std::vector<Geo::Segment2<ScalarType> > &segs)
    {
        bbox.SetEmpty();
        TMark=0;
        SegMark.clear();
        SegMark.resize(segs.size(),0);

        for (size_t i=0;i<segs.size();i++)
        {
            bbox.Add(segs[i].P(0));
            bbox.Add(segs[i].P(1));
        }

        //inflate a bit
        CoordType size3D=bbox.Max-bbox.Min;
        ScalarType inflateFact=size3D.Norm()*0.001;
        CoordType inflateVect=CoordType(inflateFact,inflateFact);
        bbox.Min=bbox.Min-inflateVect;
        bbox.Max=bbox.Max+inflateVect;
        size3D=bbox.Max-bbox.Min;

        BestDim(segs.size(),size3D,siz);

        voxel.X= size3D.X/(ScalarType)siz.X;
        voxel.Y= size3D.Y/(ScalarType)siz.Y;

        Grid.clear();

        //allocate space
        Grid.resize(siz.X);
        for (size_t i=0;i<Grid.size();i++)
        {
            Grid[i].resize(siz.Y);
        }

        for (size_t i=0;i<segs.size();i++)
            Insert(segs,i);

    }


    SegmentSpatialIndex2D()
    {TMark=0;}

};

}

#endif
