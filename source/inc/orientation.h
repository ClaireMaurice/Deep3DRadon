#pragma once

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Geometry"

#include "CImg.h"

////using namespace Eigen;
//using namespace cimg_library;

#include <iomanip>


static const double DtoR = 0.017453292519943295769236907684;
static const double RtoD = 57.295779513082320876798154814105;
static const double PI = 3.1415926535897932384626433832795;

class Euler;
class Quat;
class Rodrigues;
class AxisAnglePair;


class Miller 
{
public:
   //constructeur
   //Miller() {};
   Miller(int h=0, int k=0, int l=1, int u=1, int v=0, int w=0);
   Miller(Eigen::Vector3i hkl, Eigen::Vector3i uvw);
   Miller(const Miller &mil);

   Eigen::Vector3i Plan() const {return HKL;}
   Eigen::Vector3i Dir() const {return UVW;}

   //// operateur d'assignation
   Miller &operator = (const Miller &m);

   //// operateurs de comparaison
   //bool &operator == (const Miller &m) const;
   //bool &operator != (const Miller &m) const;

   //// fonctions de conversion
   Eigen::Matrix3d toRotationMatrix();
   // TODO : conversions indirectes
   Euler toEuler();
   // Quat toQuaternion();
   // Rodrigues toRodrigues();
   // AxisAnglePair toAxisAnglePair();

private:
    Eigen::Vector3i HKL;
    Eigen::Vector3i UVW;
};

std::ostream &operator << (std::ostream &s, const Miller &m);

class Euler
{
public:
   Euler(double phi1=0, double PHI=0, double phi2=0);
   Euler(const Euler &e);
   Euler(const Eigen::Vector3d &vec);
   Euler(const Eigen::Vector3f &vec);
   double phi1() const { return eul(0)*RtoD; }
   double phi() const { return eul(1)*RtoD; }
   double phi2() const { return eul(2)*RtoD; }
   void set_phi1(double val) { eul[0]=val; }
   void set_phi(double val) { eul[1]=val; }
   void set_phi2(double val) { eul[2]=val; }
   Eigen::Vector3d vec() const {return eul*RtoD;}

   Eigen::Matrix3d toRotationMatrix();
   Euler &fromRotationMatrix(const Eigen::Matrix3d &g);
   Euler get_fromRotationMatrix(const Eigen::Matrix3d &g);
   void asymmetric();
   Euler &smallestVariant(int spaceGroup);
private:
    Eigen::Vector3d eul;  // 3 angles d'Euler, convention de Bunge, en radians
};

std::ostream &operator << (std::ostream &s, const Euler &e);

class Quat
{
public:
    Quat(double ro=1, double lambda=0, double mu=0, double nu=0);
    Quat(double ro, Eigen::Vector3d v, bool AngleAxis=false);
    Quat(const Quat &q);
    Quat(const Euler &e);
    Quat(const Rodrigues &R);
    Quat(const AxisAnglePair &aap);
    Quat(double theta, int axisindex);

    double ro() const {return scal;}
    double lambda() const {return vec[0];}
    double mu() const {return vec[1];}
    double nu() const {return vec[2];}
    Eigen::Vector3d v() const {return vec;}
    double angle() const {return RtoD*2*acos(scal);}
    Eigen::Vector3d axis() const {return vec.normalized();}


    double &operator[](int index);
    const double &operator[](int index) const;

    // assignement
    Quat &operator = (const Quat &q);
    Quat &operator += (const Quat &q);
    Quat &operator -= (const Quat &q);
    Quat &operator *= (double alpha);
    Quat &operator *= (const Quat &q);			// attention : non commutatif
    Quat &operator /= (double alpha);

	// comparaison
    bool operator == (const Quat &q) const;
    bool operator != (const Quat &q) const;

	// operateurs arithmetiques
    Quat operator + (const Quat &q) const;
    Quat operator - (const Quat &q) const;
    Quat operator - () const;
    Quat operator * (const Quat &q) const;

	// produit scalaire et norme
	//double dot(const Quat &q1, const Quat &q2);
    Quat &Normalise();
    Quat &Conjugate();
    Quat &Inverse();
    Quat getInverse();
    Quat &sortDecreasingOrder();


	// conversion 
    Eigen::Matrix3d toRotationMatrix();
    Quat &fromRotationMatrix(const Eigen::Matrix3d &g, bool sym=false, int spacegroup=0);
    Euler toEuler();
    Quat &fromEuler(const Euler &e);
    Quat &fromRodrigues(const Rodrigues &r);
    Quat &fromAxisAnglePair(const AxisAnglePair &aap);
    Quat &fromAngleAxis(double theta, int axisindex);


    Quat &MinimalMisorientationQuat(const Quat &q1, const Quat &q2, int spacegroup=195, bool active=false);

    //CM-1908 : fonction rapide pour calculer l'angle de desorientation entre deux quaternions
    double desorientationAngle(const Quat &q2, int spaceGroup=195);

private:
    double scal;
    Eigen::Vector3d vec;
};

static const double invrac2 = 0.70710678118654752440084436210485;
static const Quat SymQuatCubic[24] = {
    Quat(1,0,0,0),
    Quat(invrac2,invrac2,0,0),
    Quat(0,1,0,0),
    Quat(-invrac2,invrac2,0,0),
    Quat(invrac2,0,invrac2,0),
    Quat(0,0,1,0),
    Quat(-invrac2,0,invrac2,0),
    Quat(invrac2,0,0,invrac2),
    Quat(0,0,0,1),
    Quat(-invrac2,0,0,invrac2),
    Quat(0,invrac2,invrac2,0),
    Quat(0,-invrac2,invrac2,0),
    Quat(0,0,invrac2,invrac2),
    Quat(0,0,-invrac2,invrac2),
    Quat(0,invrac2,0,invrac2),
    Quat(0,-invrac2,0,invrac2),
    Quat(0.5,0.5,0.5,0.5),
    Quat(0.5,-0.5,0.5,0.5),
    Quat(0.5,0.5,-0.5,0.5),
    Quat(0.5,0.5,0.5,-0.5),
    Quat(-0.5,0.5,0.5,0.5),
    Quat(-0.5,-0.5,0.5,0.5),
    Quat(-0.5,0.5,-0.5,0.5),
    Quat(-0.5,0.5,0.5,-0.5)
};

static const double halfrac3=0.86602540378443864676372317075294;
static const Quat SymQuatHexagonal[12] = {
    Quat(1,0,0,0),
    Quat(0,1,0,0),
    Quat(0,0,1,0),
    Quat(0,0,0,1),
    Quat(halfrac3,0,0,0.5),
    Quat(halfrac3,0,0,-0.5),
    Quat(0,halfrac3,0.5,0),
    Quat(0,halfrac3,-0.5,0),
    Quat(0,0.5,halfrac3,0),
    Quat(0,0.5,-halfrac3,0),
    Quat(0.5,0,0,halfrac3),
    Quat(0.5,0,0,-halfrac3)
};


std::ostream &operator << (std::ostream &s, const Quat &q);

class Rodrigues
{
public:
    Rodrigues(double R1=0.0, double R2=0.0, double R3=0.0);
    Rodrigues(const Eigen::Vector3d &v);
    Rodrigues(const Rodrigues &A);
    Rodrigues(const AxisAnglePair AAP);

    double R1() const {return R(0);}
    double R2() const {return R(1);}
    double R3() const {return R(2);}
    Eigen::Vector3d vec() const {return R;}

    Rodrigues &operator = (const Rodrigues &AAP);

    Eigen::Matrix3d toRotationMatrix();
    Rodrigues &fromAxisAnglePair();
    Rodrigues &fromQuaternion(const Quat &);
    Rodrigues &fromEigenVector(const Eigen::Vector3d &v);

private:
    Eigen::Vector3d R;
};

class AxisAnglePair
{
public:
    AxisAnglePair(double r1=1.0, double r2=0.0, double r3=0.0, double theta=0.0);
    AxisAnglePair(const Eigen::Vector3d &axis, double angle=0);
    AxisAnglePair(const AxisAnglePair &AAP);
    AxisAnglePair &operator = (const AxisAnglePair &AAP);

    double Angle() const {return RtoD * theta;}
    Eigen::Vector3d Axis() const {return r;}
    void setAngle(double val) {theta = DtoR*val;}
    void setAxis(Eigen::Vector3d &vec) {r= vec.normalized();}

    Eigen::Matrix3d toRotationMatrix();
    AxisAnglePair &fromRotationMatrix(const Eigen::Matrix3d g);

    Euler toEuler();
    AxisAnglePair &fromEuler(const Euler e);

    Rodrigues toRodrigues();
    AxisAnglePair &fromRodrigues(const Rodrigues r);

private:
    Eigen::Vector3d r;                 // vecteur unitaire de l'axe de rotation
    double theta;		// angle de rotation exprime en radians
};



#ifdef cimg_version
// la librairie peut etre utilisee sur des cartos d'angles d'Euler ou des cartos de quaternions
// on definit des macros "a la mode de CImg" pour travailler sur des voisinages

#define QImg_3x3(Q) Quat Q[9]; \
                      Quat& Q##pp = Q[0]; Quat& Q##cp = Q[1]; Quat& Q##np = Q[2]; \
                      Quat& Q##pc = Q[3]; Quat& Q##cc = Q[4]; Quat& Q##nc = Q[5]; \
                      Quat& Q##pn = Q[6]; Quat& Q##cn = Q[7]; Quat& Q##nn = Q[8]; \
                      Q##pp = Q##cp = Q##np = \
                      Q##pc = Q##cc = Q##nc = \
                      Q##pn = Q##cn = Q##nn = Quat(0,0,0,0)

#define QImg_for3x3(img,x,y,I) \
   cimg_for3((img)._height,y) for (int x = 0, \
   _p1##x = 0, \
   _n1##x = (int)( \
   (Q[0] = Q[1] = Quat((img)(_p1##x,_p1##y,0),(img)(_p1##x,_p1##y,1),(img)(_p1##x,_p1##y,2),(img)(_p1##x,_p1##y,3))), \
   (Q[3] = Q[4] = Quat((img)(0,y,0),(img)(0,y,1),(img)(0,y,2),(img)(0,y,3))), \
   (Q[6] = Q[7] = Quat((img)(0,_n1##y,0),(img)(0,_n1##y,1),(img)(0,_n1##y,2),(img)(0,_n1##y,3))),      \
   1>=(img)._width?(img).width()-1:1); \
   (_n1##x<(img).width() && ( \
   (Q[2] = Quat((img)(_n1##x,_p1##y,0),(img)(_n1##x,_p1##y,1),(img)(_n1##x,_p1##y,2),(img)(_n1##x,_p1##y,3))), \
   (Q[5] = Quat((img)(_n1##x,y,0),(img)(_n1##x,y,1),(img)(_n1##x,y,2),(img)(_n1##x,y,3))), \
   (Q[8] = Quat((img)(_n1##x,_n1##y,0),(img)(_n1##x,_n1##y,1),(img)(_n1##x,_n1##y,2),(img)(_n1##x,_n1##y,3))) ,1)) || \
   x==--_n1##x; \
   Q[0] = Q[1], Q[1] = Q[2], \
   Q[3] = Q[4], Q[4] = Q[5], \
   Q[6] = Q[7], Q[7] = Q[8], \
   _p1##x = x++, ++_n1##x)

#endif  //cimg_version