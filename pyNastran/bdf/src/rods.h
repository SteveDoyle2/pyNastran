#ifndef ROD_H
#define ROD_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

class CROD {
//    private:
    public:
        //int n1, n2;
        int eid;
        vector<int> nids;
        int pid;

    	CROD(int EID, float N1, int N2, int PID);
    	//std::string print_card_8();
    	//void print_card_16();
    	//void write();
};


class PROD {
//    private:
    public:
        //int n1, n2;
        int pid;
        int mid;
        float A;
        float J;

    	PROD(int PID, int MID, float AREA, float Ji);
    	//std::string print_card_8();
    	//void print_card_16();
    	//void write();
};

#endif
