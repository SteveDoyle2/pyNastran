#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;

#include "grid.h"

// constants
//#define PI 3.14159


//string GRID::type="GRID";

void GRID::write() {
    cout << "type  = GRID\n";
    cout << "  nid =" << nid << endl;
    cout << "  X   =" << x << endl;
    cout << "  Y   =" << y << endl;
    cout << "  Z   =" << z << endl;
    cout << "  cp  =" << cp << endl;
    cout << "  cd  =" << cd << endl;
    cout << "  seid=" << seid << endl;
}

//int GRID::getCd(void) {
//  return _cd;
//}

//GRID::GRID() {
//  _x1 = 0.0;
//  _x2 = 0.0;
//  _x3 = 0.0;
//  _cp = 0;
//  _seid = 0;
//  _cd = 0;
//}

GRID::GRID() {}
GRID::~GRID() {}

GRID::GRID(int Nid, float X, float Y, float Z, int Cp, int Cd, int Seid) {
    nid = Nid;
    cp = Cp;
    x = X;
    y = Y;
    z = Z;
    cd = Cd;
    seid = Seid;

    //_id = integer(id);
    //_x1 = float_or_blank(x1, 0.0);
    //_x2 = float_or_blank(x2, 0.0);
    //_x3 = float_or_blank(x3, 0.0);
    //_cp = integer_or_blank(cp, 0);
    //_seid = integer_or_blank(seid, 0);
    //_cd = integer_or_blank(cd, 0);
}

//string GRID::print_card_8(){
//	stringstream ss = "GRID    " << setw(8) << nid  << cp
//		<< x << y << z << cd
//		<< seid;
//	return ss.str();
//}


//Coord::GRID getCpCoord() {
//    return model.coord[cp];
//}

//int::GRID getSeid(void) {
//  return *_seid;
//}

//Coord::GRID getCdCoord()
//    return model.coord[cd];
//}


//node.getCp(model)
//node.getCp()


//GRID::~GRID () {
//  delete _x1;
//  delete _x2;
//  delete _x3;
//  delete _cp;
//  delete _seid;
//  delete _cd;
//}


/*
int main1()
{
    cout << "Hello world!\n";
    GRID g(10,1.,2.,3.,4,5,6);
    cout << "g.type="<< g.type << " g.cp=" << g.getCp() << endl;
    cout << "cd = " << g.getCd() << endl;
    //cout << "1 -> " << integer_or_blank(5) << endl;
    //cout << "2 -> " << integer_or_blank(5, 6) << endl;
    //cout << "3 -> " << integer_or_blank("5") << endl;
    //string mystring = "asdf";
    //cout << mystring << endl;
    g.write();
    //return 4;
    return 0;
};
*/

int gridmain()
{
    cout << "Hello world!\n";
    //cout << "1 -> " << integer_or_blank(5) << endl;
    //cout << "2 -> " << integer_or_blank(5, 6) << endl;
    //cout << "3 -> " << integer_or_blank("5") << endl;
    //string mystring = "asdf";
    //cout << mystring << endl;
    //g.write();
    //return 4;
    return 0;
};
