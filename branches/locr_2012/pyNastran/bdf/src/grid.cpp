#include <iostream>
#include <string>
#include <sstream>
using namespace std;

// constants
//#define PI 3.14159


int integer_or_blank(int value) {
	return value;}

//int integer_or_blank(stringstream value) {
//	return int(value);}

int integer_or_blank(int value, int default_value) {
	return default_value;}

float float_or_blank(float value, float default_value) {
	return default_value;}

class GRID {
	private:
		int _id;
		float _x1;
		float _x2;
		float _x3;
		int _cp;
		int _seid;
		int _cd;

	public:
		static string type;
		//GRID();
		GRID(int, float, float, float, int, int, int);
		int getID(void) {return _id;}
		int getCp(void) {return _cp;}
		int getSeid(void) {return _seid;}
		int getCd(void) {return _cd;}
		//~GRID();
		
		//Coord getCpCoord();
		//Coord getCdCoord();

        void write(void);
};

string GRID::type="GRID";

void GRID::write() {
	cout << "type  =" << type << endl;
	cout << "  ID  =" << _id << endl;
	cout << "  x1  =" << _x1 << endl;
	cout << "  x2  =" << _x2 << endl;
	cout << "  x3  =" << _x3 << endl;
	cout << "  cp  =" << _cp << endl;
	cout << "  cd  =" << _cd << endl;
	cout << "  seid=" << _seid << endl; 
}

//int GRID::getCd(void) {
//	return _cd;
//}

//GRID::GRID() {
//	_x1 = 0.0;
//	_x2 = 0.0;
//	_x3 = 0.0;
//	_cp = 0;
//	_seid = 0;
//	_cd = 0;
//}

GRID::GRID(int id, float x1, float x2, float x3, int cp, int seid, int cd) {
	_id = id;
	_x1 = x1;
	_x2 = x2;
	_x3 = x3;
	_cp = cp;
	_seid = seid;
	_cd = cd;

	//_id = integer(id);
	//_x1 = float_or_blank(x1, 0.0);
	//_x2 = float_or_blank(x2, 0.0);
	//_x3 = float_or_blank(x3, 0.0);
	//_cp = integer_or_blank(cp, 0);
	//_seid = integer_or_blank(seid, 0);
	//_cd = integer_or_blank(cd, 0);
}

//Coord::GRID getCpCoord() {
//    return model.coord[cp];
//}

//int::GRID getSeid(void) {
//	return *_seid;
//}

//Coord::GRID getCdCoord()
//    return model.coord[cd];
//}


//node.getCp(model)
//node.getCp()


//GRID::~GRID () { 
//	delete _x1;
//	delete _x2;
//	delete _x3;
//	delete _cp;
//	delete _seid;
//	delete _cd;
//}

int main()
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
