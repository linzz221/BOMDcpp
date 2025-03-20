#ifndef BOMD_CLASSH
#define BOMD_CLASSH

#include <string>
#include <vector>

using namespace std;
typedef vector<double> vecd;

class BOMD {
private:
	const string tempgjffile = "tempforce.gjf";
	const string templogfile = "tempforce.log";
	const string data_name = "data.txt";
	double temperature;
	string gjfname;
	string velname;
	string gjfhead;
	vector<string> atomsequ;
	vecd mass_sequ;
	vecd coor;
	vecd velocity;
	vecd force;
	long long atomnum;
	long long vsize;

	void makevel();
	void readgjf();
	void create_mass_sequ();
	void rungaussian();
	void readforce();
	void readvel();
	void create_gjf();
	void logging_data(ofstream & data, const long long cycle, const double temperature);
	double calc_temperature();

public:
	BOMD(const string& gjfname, const double& itemp = 300, 
		const string& velname = "noneed");
	~BOMD();
	void run(long long nstep, double dt = 0.5);

};

#endif // !BOMD_CLASSH
