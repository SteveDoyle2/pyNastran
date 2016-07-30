#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdlib.h>     /* atof */
using namespace std;

#include "bdf.h"
#include "grid.h"

BDF::BDF() {
	//nodes;
	//vector<Coord> coords;
}

BDF::~BDF() {}

void BDF::read_bdf(string BDF_filename) {
	bdf_filename = BDF_filename;
}
void BDF::write() {
	cout << "BDF.write()\n";
}

//---------------------------------------------------------------------------------------------

void BDF::_add_card_helper(std::vector<string> card_object, int card, string card_name) {
	if (card_name == "GRID") {
		BDF::add_grid(card_object);  // comment=comment
	}
	//else if(card_name == "CONROD") BDF::add_conrod(card_object);  // comment=comment
	//else if(card_name == "CROD") BDF::add_crod(card_object);  // comment=comment
	//else if(card_name == "PROD") BDF::add_prod(card_object);  // comment=comment
	//else if(card_name == "MAT1") BDF::add_mat1(card_object);  // comment=comment
	else {
		cout << "card_name='" << card_name << "' is not supported\n";
	}
}

int BDF::add_card(string card_lines, string card_name, string comment) {
	// return card_object
	int card_object = 1;  // wrong type
	return card_object;
}
//---------------------------------------------------------------------------------------------
/*
int integer(std::vector<string> card_object, int ifield, string name) {
	string value;
	int value_out;

	if (ifield >= card_object.size()) {
		cout << "too short...errror";
		return -9999;
	}

	value = card_object[ifield];
	return 2.0;
//	try {
//		value_out = atof(value);
//	}
//	catch (int e)
//	{
//		cout << "An exception occurred; #" << e << "\n";
//		cout << "Cannot parse '" << value << "'\n"; cout.flush();
//	}
//	return value_out;

int integer_or_blank(std::vector<string> card_object, int ifield, string name, int default_value=0) {
	string value;
	double value_out;

	if (ifield >= card_object.size()) {
		return default_value;
	}
	value = card_object[ifield];
	return 3;
//	try {
//		value_out = int(value);
//	}
//	catch (int e)
//	{
//		cout << "An exception occurred; #" << e << "\n";
//		cout << "Cannot parse '" << value << "'\n"; cout.flush();
//	}
//	return value_out;

double double_or_blank(std::vector<string> card_object, int ifield, string name, double default_value=0.0) {
	string value;
	double value_out;

	if (ifield >= card_object.size()) {
		return default_value;
	}
	value = card_object[ifield];
	return 4.0;
//	try {
//		value_out = atof(value);
//	}
//	catch (int e)
//	{
//		cout << "An exception occurred; #" << e << "\n";
//		cout << "Cannot parse '" << value << "'\n"; cout.flush();
//	}
//	return value_out;
}
*/


//---------------------------------------------------------------------------------------------

// not done...private
void BDF::add_grid(std::vector<string> card_object, string comment) {
	GRID node;
	int nid, cp;
	double x1, x2, x3;

	nid = 1;
	cp = 0;
	x1 = 1.0;
	x2 = 2.0;
	x3 = 3.0;
	//nid = integer(card_object, 1, "node_id")
	//cp = integer_or_blank(card_object, 2, "Cp")
	//x1 = double_or_blank(card_object, 3, "x1");
	//x2 = double_or_blank(card_object, 4, "x2");
	//x3 = double_or_blank(card_object, 5, "x3");

	node = GRID(nid, cp, x1, x2, x3);
	nodes.push_back(node);
}

// not done...
void BDF::add_crod(std::vector<string> card_object, string comment) {}
void BDF::add_conrod(std::vector<string> card_object, string comment) {}
void BDF::add_prod(std::vector<string> card_object, string comment) {}
void BDF::add_mat1(std::vector<string> card_object, string comment) {}
//void BDF::add_mat2(std::vector<string> card_object, string comment) {}
//void BDF::add_mat4(std::vector<string> card_object, string comment) {}
//void BDF::add_mat5(std::vector<string> card_object, string comment) {}
//void BDF::add_mat8(std::vector<string> card_object, string comment) {}
