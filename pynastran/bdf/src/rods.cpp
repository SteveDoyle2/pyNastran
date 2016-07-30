#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;
#include "rods.h"

CROD::CROD(int EID, float N1, int N2, int PID) {
	eid = EID;
	pid = PID;
	nids.push_back(N1);
	nids.push_back(N2);
}


PROD::PROD(int PID, int MID, float AREA, float Ji) {
	pid = PID;
	mid = MID;
	A = AREA;
	J = Ji;
}
