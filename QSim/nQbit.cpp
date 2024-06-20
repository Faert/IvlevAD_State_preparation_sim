#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <chrono>
#include "Qbit.h"
#include "RM_RSA.h"
#include "Alg.h"

using namespace std;

int main(int argc, char* argv[])
{
    //State preparation
    
	vector<complex<double>> vec{ 1, 2, 3, 4, 5, 6, 7, 8 };
	normalization_vec(vec); //State vec

	vector<complex<double>> res = State_preparation(vec);
	for (size_t i = 0; i < vec.size(); i++)
	{
		cout << norm(vec[i] - res[i]) << ' ';
	}
	cout << '\n';
}
