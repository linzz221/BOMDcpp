#ifndef CONSTANT_H
#define CONSTANT_H

#include <string>
#include <unordered_map>

namespace cs {
	constexpr double amu2au = 1822.888515;  //相对原子质量 -> 原子单位质量
    constexpr double temp_au2si = 3.1577464e5; //原子单位温度 -> 开尔文温度
	constexpr double coor_au2A = 0.52917721092; //原子单位长度 -> 埃
	constexpr double time_au2fs = 2.418884326505e-2; //原子单位时间 -> 飞秒
    const std::unordered_map<std::string, double> amu_mass = {
        {"H", 1.007947}, {"He", 4.0026022}, {"C", 12.01078}, {"N", 14.00672},
        {"O", 15.99943}, {"F", 18.9984032}, {"P", 30.9737612}, {"S", 32.0655},
        {"Cl", 35.4532}
        };
    }

#endif // !CONSTANT_H
