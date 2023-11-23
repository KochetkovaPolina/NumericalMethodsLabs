#include<iostream>
#include<iomanip>
using namespace std;

void gauss(double** matrixA, int N, double* dx);

double u1(double* y, double k, double a);
double u2(double* y, double k, double a);
double u3(double* y, double k, double a);

double f1(double* y_k1, double* y_k, double T_k, double k, double a) {
	return y_k1[0] - y_k[0] - T_k * u1(y_k1, k, a);
}
double f2(double* y_k1, double* y_k, double T_k, double k, double a) {
	return y_k1[1] - y_k[1] - T_k * u2(y_k1, k, a);
}
double f3(double* y_k1, double* y_k, double T_k, double k, double a) {
	return y_k1[2] - y_k[2] - T_k * u3(y_k1, k, a);
}

void Iakobian(double** J, double* y_k1, double* y_k, double T_k, double M, int n, double e1, double k, double a) {
	double dy_i = M * y_k1[0];
	double* yy = new double[n];
	if (fabs(dy_i) < e1) {
		dy_i = e1;
	}
	for (int i = 0; i < n; i++) {
		yy[i] = y_k1[i];
	}
	yy[0] = yy[0] + dy_i;
	J[0][0] = (f1(yy, y_k, T_k, k, a) - f1(y_k1, y_k, T_k, k, a)) / dy_i;
	J[1][0] = (f2(yy, y_k, T_k, k, a) - f2(y_k1, y_k, T_k, k, a)) / dy_i;
	J[2][0] = (f3(yy, y_k, T_k, k, a) - f3(y_k1, y_k, T_k, k, a)) / dy_i;


	yy[0] = yy[0] - dy_i;
	dy_i = M * y_k1[1];
	if (fabs(dy_i) < e1) {
		dy_i = e1;
	}
	yy[1] = yy[1] + dy_i;
	J[0][1] = (f1(yy, y_k, T_k, k, a) - f1(y_k1, y_k, T_k, k, a)) / dy_i;
	J[1][1] = (f2(yy, y_k, T_k, k, a) - f2(y_k1, y_k, T_k, k, a)) / dy_i;
	J[2][1] = (f3(yy, y_k, T_k, k, a) - f3(y_k1, y_k, T_k, k, a)) / dy_i;
	yy[1] = yy[1] - dy_i;
	dy_i = M * y_k1[2];
	if (fabs(dy_i) < e1) {
		dy_i = e1;
	}
	yy[2] = yy[2] + dy_i;
	J[0][2] = (f1(yy, y_k, T_k, k, a) - f1(y_k1, y_k, T_k, k, a)) / dy_i;
	J[1][2] = (f2(yy, y_k, T_k, k, a) - f2(y_k1, y_k, T_k, k, a)) / dy_i;
	J[2][2] = (f3(yy, y_k, T_k, k, a) - f3(y_k1, y_k, T_k, k, a)) / dy_i;
	delete[] yy;
}

void matrix(double** J, double** matrixA, double* y_k1, double* y_k, double T_k, double k, double a) {
	matrixA[0][0] = J[0][0];	matrixA[0][1] = J[0][1];	matrixA[0][2] = J[0][2];	matrixA[0][3] = -f1(y_k1, y_k, T_k, k, a);
	matrixA[1][0] = J[1][0];	matrixA[1][1] = J[1][1];	matrixA[1][2] = J[1][2];	matrixA[1][3] = -f2(y_k1, y_k, T_k, k, a);
	matrixA[2][0] = J[2][0];	matrixA[2][1] = J[2][1];	matrixA[2][2] = J[2][2];	matrixA[2][3] = -f3(y_k1, y_k, T_k, k, a);
}

bool Newton(double* y_k1, int n, double* Y_k, double T_k, double t_k_1, double k, double a) {
	double M = 0.01, e1 = 1e-9, e2 = 1e-9;
	double d1, d2; d1 = d2 = 0;
	int K = 0;
	int IsPrint = 0;
	int NIT = 30;
	double* yk = new double[n];
	double** J = new double* [n];
	for (int i = 0; i < n; i++) {
		J[i] = new double[n];
	}
	double* dy = new double[n];
	double** matrixA = new double* [n];
	for (int i = 0; i < n; i++) {
		matrixA[i] = new double[n + 1];
	}
	double maxF = 0, max = 0;
	double* U = new double[3];
	while (true) {
		if (IsPrint) {
			cout << setw(15) << "K=" << K << setw(15) << "d1=" << d1 << setw(15) << "d2=" << d2 << endl;
		}
		Iakobian(J, y_k1, Y_k, T_k, M, n, e1, k, a);
		matrix(J, matrixA, y_k1, Y_k, T_k, k, a);
		for (int i = 0; i < n; i++) {
			U[i] = matrixA[i][2];
		}

		gauss(matrixA, n, dy);

		for (int i = 0; i < n; i++) {
			yk[i] = y_k1[i];
		}

		for (int i = 0; i < n; i++) {
			y_k1[i] = y_k1[i] + dy[i];
		}

		if (K >= NIT) {
			return 1;
		}
		K++;

		maxF = fabs(U[0]);
		for (int i = 1; i < n; i++) {
			if (maxF < fabs(U[i])) {
				maxF = fabs(U[i]);
			}
		}
		d1 = maxF;
		for (int i = 1; i < n; i++) {
			if (fabs(y_k1[i]) < 1) {
				max = fabs(y_k1[i] - yk[i]);
			}
			else
				max = fabs((y_k1[i] - yk[i]) / yk[i]);
		}
		d2 = max;
		if (d1 <= e1 && d2 <= e2) {
			return 0;
		}

	}

	delete[] yk;
	delete[] dy;
	for (int i = 0; i < n; i++) {
		delete[] J[i];
	}
	delete[] J;
	delete[] U;
	for (int i = 0; i < n; i++) {
		delete[] matrixA[i];
	}
	delete[] matrixA;
}