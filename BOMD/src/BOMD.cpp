#include <cmath>
#include <algorithm>
#include <span>
#include <numeric>
#include <sstream>
#include <format>
#include <cstdlib>
#include <random>
#include <mkl.h>
#include <mkl_vsl.h>

#include "BOMD.h"
#include "constant.h"
#include "tools.h"

using namespace std;

BOMD::BOMD(const string & igjf, const double& itemp, const string & ivel) {
	gjfname = igjf;
	velname = ivel;
	temperature = itemp;
	readgjf();
	create_mass_sequ();

	if (velname == "noneed") {
		makevel(); 
	}
	else
	{
		readvel();
	}

}

BOMD::~BOMD() {
}

void BOMD::makevel() {
	random_device rd;
	VSLStreamStatePtr stream;
	double check = -500;
    vecd temp_v_2(vsize);
	vecd temp_v_sum(3);

	while (fabs(check - temperature) > 50){
		vslNewStream(&stream, VSL_BRNG_MT19937, rd());
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, vsize, velocity.data(), -0.5, 0.5);

		temp_v_sum = { 0, 0, 0 };
		for (auto i = 0; i < atomnum; ++i) {
            temp_v_sum[0] += velocity[3 * i];
            temp_v_sum[1] += velocity[3 * i + 1];
            temp_v_sum[2] += velocity[3 * i + 2];
		}
		cblas_dscal(3, double(1.0 / atomnum), temp_v_sum.data(), 1);

		// fs用于矫正温度
		vdMul(vsize, velocity.data(), velocity.data(), temp_v_2.data());
		auto sumv2 = cblas_ddot(vsize, mass_sequ.data(), 1, temp_v_2.data(), 1) / double(vsize);
		auto fs = sqrt(3 * (temperature / cs::temp_au2si) / sumv2);

		for (auto i = 0; i < atomnum; ++i){
			velocity[3 * i] = (velocity[3 * i] - temp_v_sum[0]) * fs;
			velocity[3 * i + 1] = (velocity[3 * i + 1] - temp_v_sum[1]) * fs;
			velocity[3 * i + 2] = (velocity[3 * i + 2] - temp_v_sum[2]) * fs;
		}
		
		check = calc_temperature();
		cout << check << endl;
	}
	vslDeleteStream(&stream);
}

void BOMD::readvel() {
	ifstream vel(velname);
    if (!vel.is_open()) {
        cerr << "Error: velocity file not found" << endl;
        exit(1);
    }
	string useless;
	double x, y, z;
	int i = 0;
	while (vel >> useless >> x >> y >> z) {
		velocity[i * 3] = x;
		velocity[i * 3 + 1] = y;
		velocity[i * 3 + 2] = z;
        ++i;
	}
	vel.close();
	if (i != atomnum) {
		cerr << "Error: velocity file is not correct" << endl;
		exit(1);
	}
}

double BOMD::calc_temperature() {
	vecd v2(vsize);
	vdMul(vsize, velocity.data(), velocity.data(), v2.data());
	double result = cblas_ddot(vsize, v2.data(), 1, mass_sequ.data(), 1);
	result = result * cs::temp_au2si / 3 / atomnum;
	return result;
}

// will init atomnum, atomsequ, coor, gjfhead, vsize
// and reserve force, velocity, mass_sequ
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
		while (iss >> element >> x >> y >> z) {
			atomsequ.push_back(element);
			coor.push_back(x / cs::coor_au2A);
			coor.push_back(y / cs::coor_au2A);
			coor.push_back(z / cs::coor_au2A);
		}		
	}
	atomnum = atomsequ.size();
	vsize = 3 * atomnum;
    force.resize(atomnum * 3);
    velocity.resize(atomnum * 3);
	mass_sequ.resize(atomnum * 3, 0);
}

// init mass_sequ
void BOMD::create_mass_sequ() {
	auto it = mass_sequ.begin();
	for (const auto& el : atomsequ) {
		double mass = cs::amu_mass.at(el) * cs::amu2au;
		fill_n(it, 3, mass);
		it += 3;
	}
}

// update force
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
	for (int i = 0; i < atomnum; ++i) {
		istringstream iss(alllogline[i]);
		iss >> useless >> useless >> x >> y >> z;
        force[i * 3] = x;
        force[i * 3 + 1] = y;
		force[i * 3 + 2] = z;
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
			gjf << format(" {:>15.10f}", coor[i * 3 + j] * cs::coor_au2A);
		}
		gjf << '\n';
	}
	gjf << "\n\n\n";
	gjf.close();
}

void BOMD::logging_data(ofstream & data ,const long long cycle, 
	const double temperature) {
	static string equ40(80, '=');
	static string sub10(10, '-');
	const string midd(format("{}    -----cycle  {} -----  Temperature {:3f}  ------ {}",
		sub10, cycle, temperature, sub10));
	data << equ40 << '\n' << midd << '\n' << equ40 << '\n';
	for (int i = 0; i < atomnum; ++i) {
		data << format("{:>4}        {:>15.10f}      {:>15.10f}      {:>15.10f}",
			atomsequ[i], coor[i * 3], coor[i * 3 + 1], coor[i * 3 + 2]) << '\n';
	}
	data << "\n";

}

void BOMD::run(const long long nstep, const double idt) {

	// 先清空data文件
	ofstream record_clear(data_name, ios::trunc);
	record_clear.close();

	ofstream record_data(data_name, ios::app);

	double dt = idt / cs::time_au2fs;
	double dt_d2 = dt / 2;
	double dt_square = dt * dt / 2;

	vecd accelerate(vsize, 0);
	vecd accelerate_old(vsize, 0);
	system(format("g16 {} {}", gjfname, templogfile).c_str());
	readforce();
	vdDiv(vsize, force.data(), mass_sequ.data(), accelerate.data());

	for (long long cycle = 1; cycle <= nstep; ++cycle) {
		// logging data
		logging_data(record_data, cycle, calc_temperature());

        // update coor
		cblas_daxpy(vsize, dt, velocity.data(), 1, coor.data(), 1);
		cblas_daxpy(vsize, dt_square, accelerate.data(), 1, coor.data(), 1);

		//update force
		create_gjf();
		rungaussian();
		readforce();
        cblas_dcopy(vsize, accelerate.data(), 1, accelerate_old.data(), 1);
		vdDiv(vsize, force.data(), mass_sequ.data(), accelerate.data());
        
		//update velocity
		cblas_daxpy(vsize, dt_d2, accelerate.data(), 1, velocity.data(), 1);
		cblas_daxpy(vsize, dt_d2, accelerate_old.data(), 1, velocity.data(), 1);

	}

}

