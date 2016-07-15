// g++ coord.cpp -o run

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstring>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::cout;
using std::endl;
using namespace std;

//#include <coord>
//#include <bdf>
#include "coord.h"
#include "grid.h"
#include "rods.h"
//#include "shell.h"
//#include "bdf.h"

int main()
{
    cout << "Hello world!\n";
    string bdfname = "fem.bdf";
//    BDF bdf(bdfname);
//    bdf.read();
//    bdf.write();
//    cout << endl;

    int cid, rid, g1, g2, g3;
    cid = 10;
    rid = 20;
    g1 = 1;
    g2 = 2;
    g3 = 3;

    CORD1R c1 = CORD1R(cid, rid, g1, g2, g3);
    CORD2R c2 = CORD2R(cid, rid,
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
