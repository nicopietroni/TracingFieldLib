#ifndef GL_DRAW_BASE
#define GL_DRAW_BASE

#include <vector>
#include <GL/glew.h>
#include "../point3.h"
#include "../point2.h"
#include <Eigen/Dense>


namespace GLDraw
{
inline void glScale(float const & p){ glScalef(p,p,p);}
inline void glScale(double const & p){ glScaled(p,p,p);}

template <class ScalarType>
inline void glColor(const Geo::Point3<ScalarType> &c)
{
    assert(c.X>=0);
    assert(c.X<=1);
    assert(c.Y>=0);
    assert(c.Y<=1);
    assert(c.Z>=0);
    assert(c.Z<=1);
    glColor3f((float)c.X,(float)c.Y,(float)c.Z);
}

template <class ScalarType>
inline void glColorAlpha(const Geo::Point3<ScalarType> &c,
                         const ScalarType &alpha)
{
    assert(c.X>=0);
    assert(c.X<=1);
    assert(c.Y>=0);
    assert(c.Y<=1);
    assert(c.Z>=0);
    assert(c.Z<=1);
    assert(alpha>=0);
    assert(alpha<=1);
    glColor4f((float)c.X,(float)c.Y,(float)c.Z,(float)alpha);
}

inline void glVertex(Geo::Point3<int> const & p)     { glVertex3i(p.X,p.Y,p.Z);}
inline void glVertex(Geo::Point3<float> const & p)   { glVertex3f(p.X,p.Y,p.Z);}
inline void glVertex(Geo::Point3<double> const & p)  { glVertex3d(p.X,p.Y,p.Z);}

inline void glVertex(Geo::Point2<int> const & p)     { glVertex3i(p.X,p.Y,0);}
inline void glVertex(Geo::Point2<float> const & p)   { glVertex3f(p.X,p.Y,0);}
inline void glVertex(Geo::Point2<double> const & p)  { glVertex3d(p.X,p.Y,0);}

inline void glNormal(Geo::Point3<float> const & p)   { glNormal3f(p.X,p.Y,p.Z);}
inline void glNormal(Geo::Point3<double> const & p)  { glNormal3d(p.X,p.Y,p.Z);}

inline void glTranslate(Geo::Point3<float> const & p) { glTranslatef(p.X,p.Y,p.Z);}
inline void glTranslate(Geo::Point3<double> const & p) { glTranslated(p.X,p.Y,p.Z);}

inline void glTranslate(Geo::Point2<float> const & p) { glTranslatef(p.X,p.Y,0);}
inline void glTranslate(Geo::Point2<double> const & p) { glTranslated(p.X,p.Y,0);}

inline void glScale(Geo::Point3<float> const & p) { glScalef(p.X,p.Y,p.Z);}
inline void glScale(Geo::Point3<double> const & p){ glScaled(p.X,p.Y,p.Z);}


template <class ScalarType>
void GLDrawSegment(const Geo::Point3<ScalarType> &p0,
                   const Geo::Point3<ScalarType> &p1,
                   const ScalarType LineWidth,
                   const ScalarType PointWidth=-1,
                   const Geo::Point3<ScalarType> colorPoint=Geo::Point3<ScalarType>(0,0,0),
                   const Geo::Point3<ScalarType> colorLine=Geo::Point3<ScalarType>(0,0,0))
{
    //    std::cout<<"P0:"<<p0.X<<","<<p0.Y<<","<<p0.Z<<std::endl;
    //    std::cout<<"P1:"<<p1.X<<","<<p1.Y<<","<<p1.Z<<std::endl;
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDepthRange(0.0,0.9999);
    glEnable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);
    if (PointWidth>0)
    {
        glColor(colorPoint);
        glPointSize (PointWidth);
        glBegin(GL_POINTS);
        glVertex (p0);
        glVertex (p1);
        glEnd ();
    }

    glColor(colorLine);
    glLineWidth (LineWidth);

    glBegin (GL_LINES);
    glVertex (p0);
    glVertex (p1);
    glEnd ();

    glPopAttrib();
}

template <class ScalarType>
void GLDrawSegment(const Geo::Point2<ScalarType> &p0,
                   const Geo::Point2<ScalarType> &p1,
                   const ScalarType LineWidth,
                   const ScalarType PointWidth=-1,
                   const Geo::Point3<ScalarType> colorPoint=Geo::Point3<ScalarType>(0,0,0),
                   const Geo::Point3<ScalarType> colorLine=Geo::Point3<ScalarType>(0,0,0))
{
    Geo::Point3<ScalarType> p0_3D(p0.X,p0.Y,0);
    Geo::Point3<ScalarType> p1_3D(p1.X,p1.Y,0);
    GLDrawSegment<ScalarType>(p0_3D,p1_3D,LineWidth,PointWidth,
                              colorPoint,colorLine);
}

template <class ScalarType>
void GLDrawPolyline(const std::vector<Geo::Point3<ScalarType> > &Poly,
                    const Geo::Point3<ScalarType> colorLine=Geo::Point3<ScalarType>(0,0,0),
                    const ScalarType LineWidth=5)
{

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDepthRange(0.0,0.9999);
    glEnable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);


    glColor(colorLine);
    glLineWidth (LineWidth);
    glBegin (GL_LINE_STRIP);

    for (size_t i=0;i<Poly.size();i++)
        glVertex (Poly[i]);

    glEnd ();

    glPopAttrib();
}

template <class ScalarType>
void GLDrawWireTriangle(const Geo::Point3<ScalarType> &P0,
                        const Geo::Point3<ScalarType> &P1,
                        const Geo::Point3<ScalarType> &P2,
                        const Geo::Point3<ScalarType> colorLine=Geo::Point3<ScalarType>(0,0,0),
                        const ScalarType LineWidth=5)
{

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDepthRange(0.0,0.9999);
    glEnable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);


    glColor(colorLine);
    glLineWidth (LineWidth);
    glBegin (GL_LINE_LOOP);

    glVertex(P0);
    glVertex(P1);
    glVertex(P2);


    glEnd ();

    glPopAttrib();
}

template <class ScalarType>
void GLDrawPolylines(std::vector<std::vector<Geo::Point3<ScalarType> > > &Polys,
                     const std::vector<Geo::Point3<ScalarType> > colorLine,
                     const ScalarType LineWidth=5)
{
    for (size_t i=0;i<Polys.size();i++)
        GLDrawPolyline(Polys[i],colorLine[i],LineWidth);
}

void GLDrawCircle (int circleStep =64)
{
    int nside = circleStep;
    const double pi2 = 3.14159265 * 2.0;
    glBegin (GL_LINE_LOOP);
    for (int i = 0; i < nside; i++) {
        glNormal3d (cos ((double)i * pi2 / nside), sin ((double)i * pi2 / nside), 0.0);
        glVertex3d (cos ((double)i * pi2 / nside), sin ((double)i * pi2 / nside), 0.0);
    }
    glEnd ();
}

template <class ScalarType>
void GLDrawOrbitSphere(const Eigen::Matrix4<ScalarType> InvMatrix,
                       const Geo::Point3<ScalarType> &center,
                       const ScalarType radius,
                       const ScalarType sca,
                       const ScalarType LineWidth)
{

    glPushAttrib(GL_TRANSFORM_BIT | GL_DEPTH_BUFFER_BIT | GL_ENABLE_BIT | GL_LINE_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix ();
    glDepthMask(GL_FALSE);

    //vcg::Point3f center = ToVCG(tb->center) + ToVCG(tb->track.InverseMatrix()*Geo::Point3<float>(0, 0, 0));
    Geo::Point3<float> centerScene = center + InvMatrix *Geo::Point3<float>(0, 0, 0);
    GLDraw::glTranslate(Geo::Point3<float>(centerScene.X,centerScene.Y,centerScene.Z));
    GLDraw::glScale (radius/sca);

    float amb[4] = { .35f, .35f, .35f, 1.0f };
    float col[4] = { .5f, .5f, .8f, 1.0f };
    glEnable (GL_LINE_SMOOTH);

    glLineWidth (LineWidth);

    glDisable(GL_COLOR_MATERIAL); // has to be disabled, it is used by wrapper to draw meshes, and prevent direct material setting, used here

    glEnable (GL_LIGHTING);
    glEnable (GL_LIGHT0);
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glColor (ColorTrack);
    glMaterialfv (GL_FRONT_AND_BACK, GL_EMISSION, amb);

    col[0] = .40f; col[1] = .40f; col[2] = .85f;
    glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, col);
    GLDraw::GLDrawCircle();

    glRotatef (90, 1, 0, 0);
    col[0] = .40f; col[1] = .85f; col[2] = .40f;
    glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, col);
    GLDraw::GLDrawCircle();

    glRotatef (90, 0, 1, 0);
    col[0] = .85f; col[1] = .40f; col[2] = .40f;
    glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, col);
    GLDraw::GLDrawCircle();

    glPopMatrix();
    glPopAttrib();
}

template <class ScalarType>
void GLDrawSphereAxis(const Eigen::Matrix4<ScalarType> InvMatrix,const Geo::Point3<ScalarType> &center,
                      const ScalarType radius,const ScalarType sca,
                      const ScalarType LineWidth)
{
    glPushAttrib(GL_TRANSFORM_BIT | GL_DEPTH_BUFFER_BIT | GL_ENABLE_BIT | GL_LINE_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glDepthMask(GL_FALSE);

    Geo::Point3<float> centerScene = center + (InvMatrix*Geo::Point3<float>(0, 0, 0));
    GLDraw::glTranslate(centerScene);
    GLDraw::glScale(radius / sca);

    float amb[4] = { .35f, .35f, .35f, 1.0f };
    float col[4] = { .5f, .5f, .8f, 1.0f };
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(LineWidth);
    glDisable(GL_COLOR_MATERIAL); // has to be disabled, it is used by wrapper to draw meshes, and prevent direct material setting, used here

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    //glColor(ColorTrack);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, amb);

    col[0] = 1.0f; col[1] = 0.0f; col[2] = 0.0f;
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col);
    glBegin(GL_LINES);
    glNormal3d(-1.0, 0.0, 0.0);
    glVertex3d(-1.2, 0.0, 0.0);
    glNormal3d( 1.0, 0.0, 0.0);
    glVertex3d( 1.2, 0.0, 0.0);
    glEnd();
    col[0] = 0.0f; col[1] = 1.0f; col[2] = 0.0f;
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col);
    glBegin(GL_LINES);
    glNormal3d(0.0,-1.0, 0.0);
    glVertex3d(0.0,-1.2, 0.0);
    glNormal3d(0.0, 1.0, 0.0);
    glVertex3d(0.0, 1.2, 0.0);
    glEnd();
    col[0] = 0.0f; col[1] = 0.0f; col[2] = 1.0f;
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, col);
    glBegin(GL_LINES);
    glNormal3d(0.0, 0.0,-1.0);
    glVertex3d(0.0, 0.0,-1.2);
    glNormal3d(0.0, 0.0, 1.0);
    glVertex3d(0.0, 0.0, 1.2);
    glEnd();

    glPopMatrix();
    glPopAttrib();
}


inline void glMultMatrix(const Eigen::Matrix4f &matrix)
{
    glMultMatrixf((const GLfloat *)(matrix.transpose().data()));
}

inline void glGetv(const GLenum  pname, Eigen::Matrix4f  & m)
{
    Eigen::Matrix4f tmp;
    glGetFloatv(pname,tmp.data());
    m = tmp;//.transpose();
}


};



#endif
