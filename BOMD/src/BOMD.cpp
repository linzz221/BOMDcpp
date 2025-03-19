#include <cmath>
#include <algorithm>
#include <span>
#include <numeric>
#include <sstream>
#include <mkl.h>
#include <mkl_vsl.h>

#include "BOMD.h"
#include "constant.h"
#include "tools.h"

using namespace std;

BOMD::BOMD(const string & igjf, const string & ivel) {
    gjfname = igjf;
    velname = ivel;
	if (velname == "noneed") { makevel(); };
	readgjf();
	create_mass_sequ();

}

void BOMD::makevel(double temperature) {
	const int size = 3 * atomnum;
	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_SFMT19937, 777);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, size, velocity.data(), 0, 0.5);
	vslDeleteStream(&stream);

	auto sumv = cblas_dasum(size, velocity.data(), 1) / size;
	auto sumv2 = cblas_ddot(size, velocity.data(), 1, velocity.data(), 1) / size;
	auto fs = sqrt(3 * temperature / sumv2);
	transform(velocity.begin(), velocity.end(), velocity.begin(), [sumv, fs](double x) {
		return (x - sumv) * fs; });
}

// will init atomnum, atomsequ, coor, gjfhead
void BOMD::readgjf() {
	vector<string> allgjfline_raw;
	readlines(gjfname, allgjfline_raw);
	auto gjf_line_num = allgjfline_raw.size();
	span<string> allgjfline(allgjfline_raw);

	for (int sl = 0; sl < gjf_line_num; ++sl) {
		if (allgjfline[sl].starts_with("#")) {
            gjfhead = accumulate(allgjfline.begin(), allgjfline.begin() + sl + 5, gjfhead, 
				[](const std::string& acc, const std::string& s) {return acc + s + '\n';});
			allgjfline = allgjfline.last(gjf_line_num - sl - 5);
			break;
		}   
	}

    for (auto el = 0; el < allgjfline.size(); ++el) {
        if (allgjfline[el].empty()) {
			allgjfline = allgjfline.first(el);
            break;
        }
    }

	string element;
	double x, y, z;
	for (const string & line : allgjfline) {
		istringstream iss(line);
        iss >> element >> x >> y >> z;
        atomsequ.push_back(element);
        coor.push_back(x);
        coor.push_back(y);
        coor.push_back(z);
	}
    atomnum = atomsequ.size();

}

// init mass_sequ
void BOMD::create_mass_sequ() {
	mass_sequ.resize(atomnum, 0);
    transform(atomsequ.begin(), atomsequ.end(), mass_sequ.begin(), [](const string& el) {
        return cs::amu_mass.at(el); });
}

void BOMD::showmass() {
	cout << gjfhead << endl;
}