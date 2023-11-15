#ifndef POINT_GRID
#define POINT_GRID

#include "grid_basics.h"
#include <mesh_interface/common_functions.h>

/*! \brief The cell of the Grid that keep the points */
template <class ScalarType>
class PointCell
{
public:
    typedef CoordType<ScalarType> CoordType;

    /*! \brief The reference to the primitives*/
    std::vector<int> entries;
    /*! \brief The position of the primitives*/
    std::vector<CoordType> positions;

    /*! \brief Base constructor */
    PointCell(){}

    /*! \brief Return the  primitives (Points) of the cell
     *  \param _entries The index of the primitive
     *  \param _positions Their position
     */
    void Get(std::vector<int> &_entries,
             std::vector<CoordType> &_positions)const
    {
        _entries=std::vector<int> (entries.begin(),entries.end());
        _positions=std::vector<CoordType> (positions.begin(),positions.end());
    }

    /*! \brief Insert a primitive in the Cell
     *  \param Elem The primitive to be inserted
     *  \param pos The position of the point
     */
    void Insert(int Elem,CoordType pos)
    {
        entries.push_back(Elem);
        positions.push_back(pos);
    }
};

/*! \brief This class implement a Grid-like search structure for points */
template <class ScalarType>
class PointSpatialIndex
{
    typedef CoordType<int> GridIndex;
    typedef PointCell<ScalarType> CellType;

    /*! \brief The actual grid of cells*/
    std::vector<std::vector<std::vector<CellType> > > Grid;

public:

    /*! \brief The bounding box of the whole space*/
   AABBox<ScalarType> bbox;
   /*! \brief Number of cells forming the grid*/
   GridIndex siz;
   /*! \brief Dimensions of a single cell*/
   CoordType<ScalarType> voxel;
   /*! \brief Global temporary Mark */


private:

   /*! \brief Return the index of the cell for a given 3D point
    *  \param p The input 3D point
    *  \param pi The grid coordinate
    */
    inline void PToIP(const CoordType<ScalarType> & p, GridIndex &pi ) const
    {
        CoordType<ScalarType> t = p - bbox.Min;
        pi.X() = int( t.X() / voxel.X() );
        pi.Y() = int( t.Y() / voxel.Y() );
        pi.Z() = int( t.Z() / voxel.Z() );
    }

    /*! \brief Convert a Bounding box to Grid Coordinates
     *  \param b The input bounding box
     *  \param ib The output Grid bounding box
     */
    inline void BoxToIBox( const AABBox<ScalarType> & b,
                           AABBox<int> & ib ) const
    {
        PToIP(b.Min, ib.Min);
        PToIP(b.Max, ib.Max);
    }

    /*! \brief Return the primitives within a cell
     *  \param Index The grid coordinates
     *  \param entries The output primitives
     *  \param positions The position of the primitives
     */
    inline void Get(const GridIndex &Index,
                    std::vector<int> &entries,
                    std::vector<CoordType<ScalarType> > &positions)const
    {
        Grid[Index.X()][Index.Y()][Index.Z()].Get(entries,positions);
    }

    /*! \brief Insert a point in a Given Cell
     *  \param Index The Grid Coordinates
     *  \param elementI The index of the element to be inserted
     *  \param Pos the coordinates of the point
     */
    void Insert(const GridIndex Index,
                const int elementI,
                CoordType<ScalarType> Pos)
    {
        Grid[Index.X()][Index.Y()][Index.Z()].Insert(elementI,Pos);
    }


    /*! \brief Insert a point the Grid
     *  \param elementI The index of the element to be inserted
     *  \param Pos the coordinates of the point
     */
    void InsertPoint(const int elementI,
                     const CoordType<ScalarType> &Pos)
    {
        //initialize temporary mark
        GridIndex IndexG;
        PToIP(Pos,IndexG);
        Insert(IndexG,elementI,Pos);
    }

public:

    /*! \brief Return the point withi a shere of a certain radius
     *  \param _p the center of the sphere
     *  \param _radius the radius of the sphere
     *  \param Elem the index of the elements inside the sphere
     *  \param Pos the positions of the elements inside the sphere
     */
    void GetInSphere(const CoordType<ScalarType> _p,
                     const ScalarType _radius,
                     std::vector<int> &Elem,
                     std::vector<CoordType<ScalarType> > &Pos)const
    {
        Elem.clear();
        Pos.clear();
        AABBox<ScalarType> RadiusBox(_p,_radius);
        AABBox<int> GridBBox;
        AABBox<int> GlobalBox(GridIndex(0,0,0),siz-GridIndex(1,1,1));
        //get the integer bbox
        BoxToIBox(RadiusBox,GridBBox);
        //intersect with global BB

        GridBBox.Intersect(GlobalBox);

        for (int ix=GridBBox.Min.X(); ix<=GridBBox.Max.X(); ix++)
            for (int iy=GridBBox.Min.Y(); iy<=GridBBox.Max.Y(); iy++)
                for (int iz=GridBBox.Min.Z(); iz<=GridBBox.Max.Z(); iz++)
                {
                    GridIndex IndexP=GridIndex(ix,iy,iz);
                    std::vector<int> test_element;
                    std::vector<CoordType<ScalarType> > positions;

                    Get(IndexP,test_element,positions);

                    for (size_t i=0;i<test_element.size();i++)
                    {
                        CoordType<ScalarType> p0=positions[i];
                        ScalarType dist=(_p-p0).Norm();
                        int elemIndex=test_element[i];
                        if (dist<_radius)
                        {
                            Elem.push_back(elemIndex);
                            Pos.push_back(p0);
                        }
                    }
                }
    }

    /*! \brief Return the closest point t oa given input position
     *  \param P0 the input position
     *  \param Clos the closest element
     *  \param MaxD the maximum distance to stop the search, if (==-1) then perform full space search
     */
    bool GetClosest(const CoordType<ScalarType> P0,
                    int &Clos,
                    ScalarType MaxD=-1)
    {
        Clos=-1;
        if (MaxD==-1)
         MaxD=(bbox.Max-bbox.Min).Norm();

        //the radius step is the size of the cell
        ScalarType step=voxel.Norm()/10;
        ScalarType curr_radius=step;
        do {
            std::vector<int> Elem;
            std::vector<CoordType<ScalarType> > Pos;
            GetInSphere(P0,curr_radius,Elem,Pos);
            ScalarType MinD=MaxD;
            if (Elem.size()!=0)
            {
                for (size_t i=0;i<Elem.size();i++)
                {
                    CoordType<ScalarType> P1=Pos[i];
                    if ((P0-P1).Norm()<MinD)
                    {
                        MinD=(P0-P1).Norm();
                        Clos=Elem[i];
                    }
                }
                return true;
            }
            curr_radius+=step;
        }while (curr_radius<MaxD);
        return false;
    }

    /*! \brief Initialize the Grid
    * \param Vertx The input vertices
    */
    void Init(const std::vector<std::vector<ScalarType> > &Vertx)
    {
        bbox.SetNull();

        for (size_t i=0;i<Vertx.size();i++)
        {
            CoordType<ScalarType> pos(Vertx[i][0],Vertx[i][1],Vertx[i][2]);
            bbox.Add(pos);
        }

        ///inflate a bit
        CoordType<ScalarType> size3D=bbox.Max-bbox.Min;
        ScalarType inflateFact=size3D.Norm()*0.001;
        CoordType<ScalarType> inflateVect=CoordType<ScalarType>(inflateFact,inflateFact,inflateFact);
        bbox.Min=bbox.Min-inflateVect;
        bbox.Max=bbox.Max+inflateVect;
        size3D=bbox.Max-bbox.Min;
        BestDim(Vertx.size(),size3D,siz);

        voxel.X()= size3D.X()/(ScalarType)siz.X();
        voxel.Y()= size3D.Y()/(ScalarType)siz.Y();
        voxel.Z()= size3D.Z()/(ScalarType)siz.Z();

        Grid.clear();

        //allocate space
        Grid.resize(siz.X());
        for (size_t i=0;i<Grid.size();i++)
        {
            Grid[i].resize(siz.Y());
            for (size_t j=0;j<Grid[i].size();j++)
                Grid[i][j].resize(siz.Z());
        }
        for (size_t i=0;i<Vertx.size();i++)
        {
            CoordType<ScalarType> pos(Vertx[i][0],Vertx[i][1],Vertx[i][2]);
            InsertPoint(i,pos);
        }
    }

};


#endif
