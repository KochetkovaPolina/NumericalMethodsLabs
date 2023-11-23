#include<iomanip>
#include<iostream>
#include<fstream>
using namespace std;

bool Newton(double* Y_k_1, int n, double* Y_k, double T_k, double t_k_1, double k, double a);

void Implicit_Euler(int n) {
	int A = 0;
	double  T = 1, t_k = 0, e_dop = 0.001, T_max = 0.1, T_min = 0.001, T_k, T_k_1, T_k1, t_k1, a = 1, k = 2, maxe_k;
	double* Y_k = new double[n];
	double* Y_k_1 = new double[n];
	double* Y_k1 = new double[n];
	double* e_k = new double[n];
	double* F = new double[n];
	for (int i = 0; i < n; i++) {
		Y_k[i] = Y_k_1[i] = Y_k1[i] = 1;
	}

	T_k_1 = T_k = T_min;
	ofstream out;
	out.open("Implicit Euler.txt", ios::trunc);
	out << "t_k" << '\t' << "Y1" << '\t' << "Y2" << '\t' << "Y3" << endl;
	while (true) {
		t_k1 = t_k + T_k;
		Newton(Y_k1, n, Y_k, T_k, t_k1, k, a);
		for (int i = 0; i < n; i++) {
			e_k[i] = -(T_k / (T_k + T_k_1)) * (Y_k1[i] - Y_k[i] - (T_k / T_k_1) * (Y_k[i] - Y_k_1[i]));
		}
		maxe_k = fabs(e_k[0]);
		for (int i = 1; i < n; i++) {
			if (fabs(e_k[i]) > maxe_k)
				maxe_k = fabs(e_k[i]);
		}
		if (maxe_k > e_dop) {
			T_k = T_k / 2;
			t_k1 = t_k;
			for (int i = 0; i < n; i++) {
				Y_k1[i] = Y_k[i];
			}
			A++;
			continue;
		}	
		T_k1 = pow((e_dop / maxe_k), 0.5) * T_k;
		if (T_k1 > T_max) {
			T_k1 = T_max;
		}

		
		/*if (maxe_k <= e_dop && e_dop / 4 < maxe_k) {
			T_k1 = T_k;
		}
		if (maxe_k <= e_dop / 4) {
			T_k1 = 2 * T_k;
		}

		if (T_k1 > T_max) {
			T_k1 = T_max;
		}*/

		out << t_k1 << '\t' << Y_k1[0] << '\t' << Y_k1[1] << '\t' << Y_k1[2] << endl;
		for (int i = 0; i < n; i++) {
			Y_k_1[i] = Y_k[i];
			Y_k[i] = Y_k1[i];
		}
		T_k_1 = T_k;
		T_k = T_k1;
		t_k = t_k1;
		if (t_k >= T)
			break;
	}
	//cout << "A=" << A << endl;
	out.close();
	delete[] Y_k;
	delete[] Y_k_1;
	delete[] Y_k1;
	delete[] F;
	delete[] e_k;
}