#ifndef _COORD_

#define _COORD_

// g++ coord.cpp -o run

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

//#include <coord>
//#include <bdf>
#include "coord.h"
#include "bdf.h"

// constants
//#define PI 3.14159

/*
array normalize(array v) {
    float normV = sqrt(v1**2 + v2**2 + v3**2);
    return v/normV;
}
*/

//------------------------------------------------------------------------
Coord::Coord(int CID, int RID) {
    cid = CID;
    rid = RID;
    //cid = integer(CID);
    //rid = integer_or_blank(RID);
}

/*
void Coord::setup() {
    e13 = e3 - e1;
    e12 = e2 - e1;
    k = normalize(e12);
    j = normalize(cross(k, e13));
    i = cross(j, k);
}

Coord::transform_to_local(array p, Matrix matrix) {
    pCoord = dot(p - e1, transpose(matrix));
    pLocal = XYZtoCoord(pCoord);
    return pLocal;
}
*/
//void Coord::write() {
//  cout << "Coord write method not implemented\n";
//}

//------------------------------------------------------------------------
class RectangularCoord {
    //public:
        //RectangularCoord(int, int);
        //XYZtoCoord(array);
        //coordToXYZ(array);
};

/*
array RectangularCoord::XYZtoCoord(array p) {
    return p; }

array RectangularCoord::coordToXYZ(array p) {
    return p; }
*/

//------------------------------------------------------------------------
class CylindricalCoord {
    //public:
        //CylindricalCoord(int, int);
        //XYZtoCoord(array);
        //coordToXYZ(array);
};

/*
array CylindricalCoord::XYZtoCoord(array p) {
    R = p[0];
    theta = radians(p[1]);
    x = R * cos(theta);
    y = R * sin(theta);
    return array(x, y, p[2]) + e1;
}

array CylindricalCoord::coordToXYZ(array p) {
    (x, y, z) = p;
    theta = degrees(atan2(y, x));
    R = sqrt(x*x + y*y);
    return array([R, theta, z])
}
*/

//------------------------------------------------------------------------
class SphericalCoord {
    public:
        //SphericalCoord(int, int);
        //XYZtoCoord(array);
        //coordToXYZ(array);
};

/*
array SphericalCoord::XYZtoCoord(array p) {
    (x, y, z) = p;
    R = sqrt(x*x + y*y + z*z);
    theta = degrees(atan2(y,x);
    if R > 0:
        phi = degrees(acos(z / R));
    else:
        phi = 0.0;
    return array([R, theta, phi])
}

array SphericalCoord::coordToXYZ(array p) {
    R = p[0];
    theta = radians(p[1]);
    phi = radians(p[2]);
    x = R * cos(theta) * sin(phi);
    y = R * sin(theta) * cos(phi);
    z = R * cos(phi);
    return array([x, y, z]) + e1
}

*/

//------------------------------------------------------------------------
class COORD1R:  public Coord {
    private:
        //int cid;
        //int rid;
        float g1;
        float g2;
        float g3;

    public:
        static string type;
        //COORD();
        COORD1R(int cid, int rid,
                int G1, int G2, int G3);
        int getG1(void) {return g1;}
        int getG2(void) {return g2;}
        int getG3(void) {return g3;}

        Coord getLocation1(BDF model, int cid=0);
        Coord getLocation2(BDF model, int cid=0);
        Coord getLocation3(BDF model, int cid=0);
        void write(void);
        //~COORD1R();

        //Coord getRidCoord();
};

string COORD1R::type="COORD1R";

COORD1R::COORD1R(int CID, int RID,
                 int G1,  int G2, int G3) : Coord(CID, RID) {
    g1 = G1;
    g2 = G2;
    g3 = G3;
}

/*
Coord COORD1R::getLocation1(BDF model, int cid) {
    return model.coords[g1].Position(cid); }

Coord COORD1R::getLocation2(BDF model, int cid) {
    return model.coords[g2].Position(cid); }

Coord COORD1R::getLocation3(BDF model, int cid) {
    return model.coords[g3].Position(cid); }
*/

void COORD1R::write() {
    cout << "type  =" << type << endl;
    cout << "  cid =" << cid << endl;
    cout << "  rid =" << rid << endl;
    cout << "  g1  =" << g1  << endl;
    cout << "  g2  =" << g2  << endl;
    cout << "  g3  =" << g3  << endl;
}

//COORD1R::~COORD1R () {
//  delete cid;
//  delete rid;
//  delete g1;
//  delete g2;
//  delete g3;
//}

//------------------------------------------------------------------------

class COORD2R  : public Coord {
    private:
        float _a1;
        float _a2;
        float _a3;

        float _b1;
        float _b2;
        float _b3;

        float _c1;
        float _c2;
        float _c3;

    public:
        static string type;
        COORD2R(int cid, int rid,
                float a1, float a2, float a3,
                float b1, float b2, float b3,
                float c1, float c2, float c3);
        int getA(void) {return _a1;}
        int getB(void) {return _b1;}
        int getC(void) {return _c1;}

        //Coord getRidCoord();
        void write(void);
        //~COORD2R();
};

string COORD2R::type="COORD2R";

void COORD2R::write() {
    cout << "type  =" << type << endl;
    cout << "  cid =" << cid << endl;
    cout << "  rid =" << rid << endl;

    cout << "  a1  =" << _a1 << endl;
    cout << "  a2  =" << _a2 << endl;
    cout << "  a3  =" << _a3 << endl;

    cout << "  b1  =" << _b1 << endl;
    cout << "  b2  =" << _b2 << endl;
    cout << "  b3  =" << _b3 << endl;

    cout << "  c1  =" << _c1 << endl;
    cout << "  c2  =" << _c2 << endl;
    cout << "  c3  =" << _c3 << endl;
}

//int COORD2R::getCd(void) {
//  return _cd;
//}

COORD2R::COORD2R(int CID, int RID,
                float a1, float a2, float a3,
                float b1, float b2, float b3,
                float c1, float c2, float c3) : Coord(CID, RID) {
    _a1 = a1;  _a2 = a2;  _a3 = a3;
    _b1 = b1;  _b2 = b2;  _b3 = b3;
    _c1 = c1;  _c2 = c2;  _c3 = c3;
}

//Coord::COORD getCpCoord() {
//    return model.coord[cp];
//}

//node.getCp(model)
//node.getCp()


//COORD2R::~COORD2R () {
//  delete cid;
//  delete rid;
//  delete _a1;  delete a2;  delete a3;
//  delete _b1;  delete b2;  delete b3;
//  delete _c1;  delete c2;  delete c3;
//}

int main()
{
    cout << "Hello world!\n";
    string bdfname = "fem.bdf";
    BDF bdf(bdfname);
    bdf.read();
    bdf.write();
    cout << endl;

    int cid, rid, g1, g2, g3;
    cid = 10;
    rid = 20;
    g1 = 1;
    g2 = 2;
    g3 = 3;

    COORD1R c1(cid, rid, g1, g2, g3);
    COORD2R c2(cid, rid,
               1.1, 2.2, 3.3,
               4.4, 5.5, 6.6,
               7.7, 8.8, 9.9);
    //cout << "g.type="<< g.type << " g.cp=" << g.getCp() << endl;
    //cout << "cd = " << g.getCd() << endl;
    c1.write();
    cout << endl;
    c2.write();
    //cout << "1 -> " << integer_or_blank(5) << endl;
    //cout << "2 -> " << integer_or_blank(5, 6) << endl;
    //cout << "3 -> " << integer_or_blank("5") << endl;
    //string mystring = "asdf";
    //cout << mystring << endl;
    //g.write();
    //return 4;
    return 0;
};
