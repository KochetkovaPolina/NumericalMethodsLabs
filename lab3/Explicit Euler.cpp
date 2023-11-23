#include<iostream>
#include<iomanip>
#include<fstream>
using namespace std;


double u1(double* u, double k, double a) {  //производные
	return (((k - a) / a) * u[1] * u[2]);
}
double u2(double* u, double k, double a) {
	return  (((k + a) / k) * u[0] * u[2]);
}
double u3(double* u, double k, double a) {
	return  (((a - k) / a) * u[0] * u[1]);
}
void VectorF(double* F, double* u, double k, double a) {
	F[0] = u1(u, k, a);
	F[1] = u2(u, k, a);
	F[2] = u3(u, k, a);
}


void Explicit_Euler(int n) {
	double  T = 1, t_k = 0, e_dop = 0.001, T_max = 0.01, T_min = 0, a = 1, k = 2;
	double* T_k = new double[n];
	double* Y_k = new double[n];
	double* F = new double[n];
	for (int i = 0; i < n; i++) {
		Y_k[i] = 1;
	}
	ofstream out;
	out.open("Explicit Euler.txt", ios::trunc);
	out << "t_k" << '\t' << "Y1" << '\t' << "Y2" << '\t' << "Y3" << endl;
	while (t_k < T) {
		VectorF(F, Y_k, k, a);
		for (int i = 0; i < n; i++) {
			T_k[i] = e_dop / (fabs(F[i]) + (e_dop / T_max)); //(3.11)
		}
		T_min = T_k[0];
		for (int i = 1; i < n; i++) {
			if (T_min > T_k[i])
				T_min = T_k[i];
		}
		for (int i = 0; i < n; i++) {
			Y_k[i] += T_min * F[i];
		}
		t_k += T_min;
		out << t_k;
		for (int i = 0; i < n; i++) {
			out << '\t' << Y_k[i];
		}
		out << endl;
	}
	out.close();
	delete[] F;
	delete[] Y_k;
	delete[] T_k;
}