#ifndef BDF_H
#define BDF_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

#include "grid.h"
//#include "coords.h"

class BDF {
    private:
    	void _add_card_helper(std::vector<string> card_object, int card, string card_name);

		// not done...
		void add_grid(std::vector<string> card_object, string comment="");
		void add_crod(std::vector<string> card_object, string comment="");
		void add_conrod(std::vector<string> card_object, string comment="");
		void add_prod(std::vector<string> card_object, string comment="");
		void add_mat1(std::vector<string> card_object, string comment="");

    public:
    	string bdf_filename;
    	std::vector<GRID> nodes;
    	//std::vector<Coord> coords;
    	//std::vector<Property> properties;
    	//std::vector<Material> materials;

    	BDF();
    	~BDF();
    	void read_bdf(string BDF_filename);
    	int  add_card(string card_lines, string card_name, string comment="");  // wrong type
    	void write();

};


#endif
