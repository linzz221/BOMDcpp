#include <cmath>
#include <algorithm>
#include <span>
#include <numeric>
#include <sstream>
#include <format>
#include <cstdlib>
#include <mkl.h>
#include <mkl_vsl.h>

#include "BOMD.h"
#include "constant.h"
#include "tools.h"

using namespace std;

BOMD::BOMD(const string & igjf, const string & ivel) {
    gjfname = igjf;
    velname = ivel;
	readgjf();
	create_mass_sequ();
	if (velname == "noneed") { makevel(); };

}

void BOMD::makevel(const double& temperature) {
	const int size = 3 * atomnum;
	velocity.resize(size, 0);
	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_MT19937, 520);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, size, velocity.data(), 0, 0.5);
	vslDeleteStream(&stream);

	auto sumv = cblas_dasum(size, velocity.data(), 1) / size;
	vdMul(size, velocity.data(), velocity.data(), velocity.data());
	auto sumv2 = cblas_ddot(size, mass_sequ.data(), 1, velocity.data(), 1) / size;
	auto fs = sqrt(3 * (temperature / cs::temp_au2si) / sumv2);
	transform(velocity.begin(), velocity.end(), velocity.begin(), [sumv, fs](double x) {return (x - sumv) * fs; });
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
        coor.push_back(x / cs::coor_au2A);
        coor.push_back(y / cs::coor_au2A);
        coor.push_back(z / cs::coor_au2A);
	}
    atomnum = atomsequ.size();

}

// init mass_sequ
void BOMD::create_mass_sequ() {
    mass_sequ.resize(atomnum * 3, 0);
    auto it = mass_sequ.begin();
    for (const auto& el : atomsequ) {
        double mass = cs::amu_mass.at(el) * cs::amu2au;
        fill_n(it, 3, mass);
        it += 3;
    }
}

void BOMD::readforce() {
	vector<string> alllogline_raw;
	const string findflag = " Center     Atomic                   Forces (Hartrees/Bohr)";
	readlines(templogfile, alllogline_raw);

	span<string> alllogline(alllogline_raw);
	int linesize = alllogline.size();
	for (int cl = 0;cl < linesize;++cl) {
		if (alllogline[cl].starts_with(findflag)) {
			alllogline = alllogline.subspan(cl+3, atomnum);
			break;
		}
	}
	
	int useless;
	double x, y, z;
	for (const string& i : alllogline) {
		istringstream iss(i);
		iss >> useless >> useless >> x >> y >> z;
		force.push_back(x);
		force.push_back(y);
		force.push_back(z);
	}

}

void BOMD::rungaussian() {
	system(format("g16 {} {}", tempgjffile, templogfile).c_str());
}

void BOMD::create_gjf() {
	ofstream gjf(tempgjffile);
	gjf << gjfhead;
	for (int i = 0; i < atomnum; ++i) {
		gjf << atomsequ[i];
		for (int j = 0; j < 3; ++j) {
			gjf << format(" {:>15.10f}", coor[i * 3 + j]);
		}
		gjf << '\n';
	}
	gjf << "\n\n\n";
	gjf.close();
}

void BOMD::run(const long long nstep, const double idt) {
	double dt = idt / cs::temp_au2si;
	double dt2 = dt * dt / 2;
}

void BOMD::showmass() {
	create_gjf();
}