#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

class GRID {
//    private:
    public:
        int nid, cp, cd;
        float x, y, z;
        int seid;
        int spc;


    	GRID(int id, float x1, float x2, float x3, int cp=0, int cd=0, int seid=0);
    	//std::string print_card_8();
    	//void print_card_16();
    	void write();
};


#endif
