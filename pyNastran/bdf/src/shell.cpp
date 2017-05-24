#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
using namespace std;
#include "shell.h"

ShellElement::ShellElement(int EID, int PID) {
	eid = EID;
	pid = PID;
}

CTRIA3::CTRIA3(int EID, int PID, int N1, int N2, int N3, int MCID) {
	eid = EID;
	pid = PID;
	nids.push_back(N1);
	nids.push_back(N2);
	nids.push_back(N3);
	mcid = MCID;
	theta = 0.;
	is_mcid = 1;
}

CTRIA3::CTRIA3(int EID, int PID, int N1, int N2, int N3, float THETA) {
	eid = EID;
	pid = PID;
	std::vector<int> nids;
	nids.push_back(N1);
	nids.push_back(N2);
	nids.push_back(N3);
	mcid = 0;
	theta = THETA;
	is_mcid = 0;
}

CQUAD4::CQUAD4(int EID, int PID, int N1, int N2, int N3, int N4, int MCID) {
	eid = EID;
	pid = PID;
	nids.push_back(N1);
	nids.push_back(N2);
	nids.push_back(N3);
	nids.push_back(N4);
	mcid = MCID;
	theta = 0.;
	is_mcid = 1;
}

CQUAD4::CQUAD4(int EID, int PID, int N1, int N2, int N3, int N4, float THETA) {
	eid = EID;
	pid = PID;
	std::vector<int> nids;
	nids.push_back(N1);
	nids.push_back(N2);
	nids.push_back(N3);
	nids.push_back(N4);
	mcid = 0;
	theta = THETA;
	is_mcid = 0;
}

PSHELL::PSHELL(int PID, int MID1, float T, int MID2, float BMIR, int MID3, float TST,
		   float Z1, float Z2, int MID4) {
	pid = PID;
	mid1 = MID1;
	mid2 = MID2;
	mid3 = MID3;
	mid4 = MID4;
	z1 = Z1;
	z2 = Z2;
	bmir = BMIR;
	tst = TST;
	nsm = 0.0;
}

/*
PCOMP::PCOMP(int PID, std::vector<int> MIDs, std::vector<float> THETAs, std::vector<int> SOUTs) {
	pid = PID;
	mids = MIDs;
	thetas = THETAs;
	souts = SOUTs;
	z1 = 0.0;
	z2 = 0.0;
	nsm = 0.0;
}
*/
