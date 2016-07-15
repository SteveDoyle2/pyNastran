#ifndef COORD_H
#define COORD_H

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

class Coord {
//private:
    //Type;
//public:
//	void write();
};

//------------------------------------------------------------------------

class CylindricalCoord {
    //public:
        //CylindricalCoord(int, int);
        //XYZtoCoord(array);
        //coordToXYZ(array);
};

class RectangularCoord {
    //public:
        //RectangularCoord(int, int);
        //XYZtoCoord(array);
        //coordToXYZ(array);
};

class SphericalCoord {
    public:
        //SphericalCoord(int, int);
        //XYZtoCoord(array);
        //coordToXYZ(array);
};

//------------------------------------------------------------------------

class CORD1R:  public Coord {
    private:
        int cid;
        int rid;
        float g1;
        float g2;
        float g3;

    public:
        static string type;
        //COORD();
        CORD1R(int cid, int rid,
               int G1, int G2, int G3);
        //int getG1(void) {return g1;}
        //int getG2(void) {return g2;}
        //int getG3(void) {return g3;}

        //Coord getLocation1(BDF model, int cid=0);
        //Coord getLocation2(BDF model, int cid=0);
        //Coord getLocation3(BDF model, int cid=0);

        void write();
        //~CORD1R();

        //Coord getRidCoord();
};

//------------------------------------------------------------------------

class CORD2R  : public Coord {
    private:
        int cid;
        int rid;

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
        CORD2R(int cid, int rid,
               float a1, float a2, float a3,
               float b1, float b2, float b3,
               float c1, float c2, float c3);
        int getA(void) {return _a1;}
        int getB(void) {return _b1;}
        int getC(void) {return _c1;}

        //Coord getRidCoord();
        void write();
        //~CORD2R();
};

#endif
