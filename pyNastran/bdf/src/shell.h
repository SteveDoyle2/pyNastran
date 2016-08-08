#ifndef SHELL_H
#define SHELL_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

class ShellElement {
	public:
		int eid, pid;
		ShellElement(int EID, int PID);
};

// : ShellElement(eid, pid)
class CTRIA3 {
//    private:
    public:
        //int n1, n2;
        int eid, pid;
        std::vector<int> nids;
        int is_mcid;
        int mcid;
        float theta;

    	CTRIA3(int EID, int PID, int N1, int N2, int N3, float THETA=0.0);
    	CTRIA3(int EID, int PID, int N1, int N2, int N3, int MCID);
    	//std::string print_card_8();
    	//void print_card_16();
    	//void write();
};



class CQUAD4 {
//    private:
    public:
        //int n1, n2;
        int eid, pid;
        std::vector<int> nids;
        int is_mcid;
        int mcid;
        float theta;

    	CQUAD4(int EID, int PID, int N1, int N2, int N3, int N4, float THETA=0.0);
    	CQUAD4(int EID, int PID, int N1, int N2, int N3, int N4, int MCID);
    	//std::string print_card_8();
    	//void print_card_16();
    	//void write();
};

class PSHELL {
	public:
		int pid, mid1, mid2, mid3, mid4;
		float t, tst, bmir;
		float nsm, z1, z2;

		PSHELL(int PID, int MID1, float T, int MID2=-1, float BMIR=1.0, int MID3=0, float TST=0.833333, float Z1=-1.0, float Z2=-1.0, int MID4=0);
		//void write();
};

/*
class PCOMP {
	public:
		int pid;
		std::vector<int> mids;
		std::vector<int> souts;
		std::vector<float> thetas;
		float nsm, z1, z2;

		PCOMP(int PID, std::vector<int> MIDs, std::vector<float> THETAs, std::vector<int> SOUTs) {
		//void write();
};
*/


#endif
