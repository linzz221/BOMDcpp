#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

void readlines(const string& filename, vector<string>& result) {
	ifstream file(filename);
	if (!file.is_open()) {
		cerr << "Error: file not found" << endl;
		exit(1);
	}
	string line;
	while (getline(file, line)) {
		result.push_back(line);
	}
	file.close();
};


template <typename T> void showvector(const vector<T>& vec) {
    for (const auto& el : vec) {
        cout << el << endl;
    }
	cout << endl;
};


#endif // !TOOLS_H
