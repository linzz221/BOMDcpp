#ifndef BOMD_CLASSH
#define BOMD_CLASSH

#include <string>
#include <vector>

using namespace std;
typedef vector<double> vecd;

class BOMD {
private:
	string gjfname;
	string velname;
	string gjfhead;
	vector<string> atomsequ;
	vecd mass_sequ;
	vecd coor;
	vecd velocity;
	vecd force;
	int atomnum;

	void makevel(double temperature = 300);
    void readgjf();
	void create_mass_sequ();
	void rungaussian();
	void readforce();
	void readvel();
	void create_gjf();
	void logging_data();

public:
	BOMD(const string& gjfname, const string& velname = "noneed");
	void showmass();
	void run(int nstep, double dt = 0.5);

};

#endif // !BOMD_CLASSH
