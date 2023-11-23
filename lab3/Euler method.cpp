//#include<iostream>
//#include<iomanip>
//#include"Euler method.h"
//
//using namespace std;
//
//const int N = 3; const double T = 1, TAU_MIN = 0.001, TAU_MAX = 0.1, LOCAL_ERROR[2] = { 1e-3,1e-5 };
//
//void Newton_method(double** x, double* dx, double** a, double* b, double tau, double** jacobian, double M, double n,
//	double (foo1)(double**, double*, double*, int), double (foo2)(double**, double*, double*, int), double (foo3)(double**, double*, double*, int),
//	double e1 = 1e-9, double e2 = 1e-9);
//
//void main() {
//	double*** h_coef_a = new double** [N];
//	double h[N] = { -6,-12,-18 };
//
//	for (int i = 0; i < N; i++) {
//		h_coef_a[i] = new double* [N];
//		for (int j = 0; j < N; j++) {
//			h_coef_a[i][j] = new double[N];
//		}
//	}
//	double** h_coef_b = new double* [N];
//	for (int i = 0; i < N; i++) {
//		h_coef_b[i] = new double[N];
//	}
//	double** a = new double* [N];
//	for (int i = 0; i < N; i++) {
//		a[i] = new double[N];
//	}
//	double* b = new double[N];
//	for (int i = 0; i < N; i++) {
//		b[i] = 0;
//		for (int j = 0; j < N; j++) {
//			a[i][j] = 0;
//		}
//	}
//	double** u = new double* [N];
//	for (int i = 0; i < N; i++) {
//		u[i] = new double[N];
//	}
//
//
//	ifstream fin("h_coef.txt");
//	readFromFile(h_coef_a, fin, N);
//	readFromFile(h_coef_b, fin, N);
//	fin.close();
//
//	filling(a, b, h_coef_a, h_coef_b, h, N);
//
//	implictEuler(u, a, b, foo1, foo2, foo3, TAU_MAX, TAU_MIN, LOCAL_ERROR[0], T, N);
//	u[0][0] = 10; u[0][1] = 22; u[0][2] = 9;
//	explictEuler(u[0], a, b, foo1, foo2, foo3, LOCAL_ERROR[0], TAU_MAX, T, N);
//
//
//	for (int i = 0; i < N; i++) {
//		for (int j = 0; j < N; j++) {
//			delete[] h_coef_a[i][j];
//		}
//		delete[] h_coef_a[i];
//	}
//	delete[] h_coef_a;
//	h_coef_a = nullptr;
//	for (int i = 0; i < N; i++) {
//		delete[] h_coef_b[i];
//	}
//	delete[] h_coef_b;
//	h_coef_b = nullptr;
//	for (int i = 0; i < N; i++) {
//		delete[] a[i];
//	}
//	delete[] a;
//	a = nullptr;
//	delete[] b;
//	b = nullptr;
//	for (int i = 0; i < N; i++) {
//		delete[] u[i];
//	}
//	delete[] u;
//	u = nullptr;
//	system("pause");
//}
//
//
//double foo1(double** a, double* b, double* u, int n) {
//	double result = 0;
//	for (int i = 0; i < n; i++) {
//		result += a[0][i] * u[i];
//	}
//	return result - b[0];
//}
//
//double foo2(double** a, double* b, double* u, int n) {
//	double result = 0;
//	for (int i = 0; i < n; i++) {
//		result += a[1][i] * u[i];
//	}
//	return result - b[1];
//}
//
//double foo3(double** a, double* b, double* u, int n) {
//	double result = 0;
//	for (int i = 0; i < n; i++) {
//		result += a[2][i] * u[i];
//	}
//	return result - b[2];
//}
//
//void readFromFile(double*** a, istream& fin, int n) {
//	double temp = 0;
//	for (int i = 0; i < n; i++) {// чтение из файла в матрицу(систему)
//		for (int j = 0; j < n; j++) {
//			for (int k = 0; k < n; k++) {
//				fin >> temp;
//				a[i][j][k] = temp;
//			}
//		}
//	}
//}
//
//void readFromFile(double** a, istream& fin, int n) {
//	double temp = 0;
//	for (int i = 0; i < n; i++) {// чтение из файла в матрицу(систему)
//		for (int j = 0; j < n; j++) {
//			fin >> temp;
//			a[i][j] = temp;
//		}
//	}
//}
//
//void filling(double** a, double* b, double*** h_coef_a, double** h_coef_b, double* h, int n) {
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < n; j++) {
//			for (int k = 0; k < n; k++) {
//				a[i][j] += (h[k] * h_coef_a[i][j][k]) / 6;
//			}
//			b[i] += (h[j] * h_coef_b[i][j]) / (6);
//		}
//	}
//
//}
//
//template<typename T>
//void print(T*** a, int n) {
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < n; j++) {
//			for (int k = 0; k < n; k++) {
//				cout << setw(10) << a[i][j][k] << " ";
//			}
//			cout << "	";
//		}
//		cout << endl;
//	}
//}
//
//template<typename T>
//void print(T** a, int n) {
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < n; j++) {
//			cout << setw(10) << a[i][j] << " ";
//		}
//		cout << "\n";
//	}
//}
//
//template<typename T>
//void print(T* a, int n) {
//	for (int i = 0; i < n; i++) {
//		cout << "x" << i << "= " << setw(6) << a[i] << "\n";
//	}
//}
//
//void implictEuler(double** y, double** a, double* b, double (*foo1) (double**, double*, double*, int),
//	double (*foo2) (double**, double*, double*, int),
//	double (*foo3) (double**, double*, double*, int),
//	double tau_max, double tau_min, double local_error, double T, int n) {
//	double* t = new double[2];	//t[0]=t[k], t[1]=t[k+1]
//	double* tau = new double[n];//tau[0]=tau[k-1], tau[1]=tau[k], tau[2]=tau[k+1]
//	double* dx = new double[n];
//	double** jacobian = new double* [n];
//	for (int i = 0; i < n; i++) {
//		jacobian[i] = new double[n];
//	}
//	double* eps = new double[n];
//	for (int i = 0; i < n; i++) {
//		eps[i] = 0;
//	}
//	typedef double (*funcmas_t) (double** a, double* b, double* u, int n);
//	funcmas_t  functions[3];
//	functions[0] = foo1;
//	functions[1] = foo2;
//	functions[2] = foo3;
//	t[0] = 0;
//	for (int i = 0; i < 3; i++) {		//пункт 2
//		y[i][0] = 10;
//		y[i][1] = 22;
//		y[i][2] = 9;
//		tau[i] = tau_min;
//	}
//	int it_num = 0;
//	double min = 0;
//	ofstream fout("Euler method.txt");
//	fout << "implict Euler:\nt	u1	u2	u3:\n";
//	while (t[1] < T) {
//
//
//		t[1] = t[0] + tau[1];				//пункт 3	-	шаг
//
//		Newton_method(y, dx, a, b, tau[1], jacobian, tau[1], n, foo1, foo2, foo3);	//пункт 4	-	решение F(y[k+1])=y[k+1]-y[k]-tau[k]*f[k](y[k+1],t[k+1])
//
//		/*cout << "\nU[k+1]:";
//		print(y[2], n);*/
//		for (int i = 0; i < n; i++) {
//			eps[i] = -(tau[1] / (tau[1] + tau[0])) * (y[2][i] - y[1][i] - (tau[1] / tau[0]) * (y[1][i] - y[0][i]));//пункт 5	-	Вычисление eps
//
//		}
//		min = eps[0];
//		for (int i = 1; i < n; i++) {
//			if (min > eps[i]) {
//				min = eps[i];
//			}
//		}
//		if (fabs(min) > local_error) {
//			tau[1] /= 2;
//			t[1] = t[0];
//			for (int i = 0; i < n; i++) {
//				y[2][i] = y[1][i];
//			}
//			continue;
//		}				//пункт 6
//
//		tau[2] = (fabs(eps[0]) > local_error) ? tau[1] / 2 : (fabs(eps[0]) > local_error / 4 && fabs(eps[0]) <= local_error) ? tau[1] : tau[1] * 2;
//		for (int i = 1; i < n; i++) {
//			min = (fabs(eps[i]) > local_error) ? tau[1] / 2 : (fabs(eps[i]) > local_error / 4 && fabs(eps[i]) <= local_error) ? tau[1] : tau[1] * 2;
//			if (min < tau[2]) tau[2] = min;
//		}			//пункт 7 ПРАВИЛО 3-Х ЗОН
//
//		/*tau[2] = sqrtf(local_error / eps[0]) * tau[1];
//		for (int i = 1; i < n; i++) {
//			min = sqrtf(local_error / eps[i]) * tau[1];
//			if (min < tau[2]) tau[2] = min;
//		}			//пункт 7 КВАЗИОПТИМАЛЬНЫЙ ВЫБОР ШАГА	*/
//
//		if (tau[2] > tau_max) tau[2] = tau_max; //пункт 8
//		/*cout << "\n\ny[" << it_num << "]\n";
//		print(y[2], n);
//		cout << "\n\nt[" << it_num << "]\n"<<t[1];		//пункт 9 -	Вывод на печать*/
//		it_num++;
//
//		fout << t[1];
//		for (int i = 0; i < n; i++) {
//			fout << "	" << y[2][i];
//		}
//		fout << "\n";//запись в файл
//
//		for (int i = 0; i < n; i++) {
//			y[0][i] = y[1][i]; y[1][i] = y[2][i];
//		}
//		tau[0] = tau[1]; tau[1] = tau[2]; t[0] = t[1];//пункт 10 - Сдвиг переменных
//	}
//	fout.close();
//
//
//
//
//	delete[] t;
//	t = nullptr;
//	delete[] tau;
//	tau = nullptr;
//	delete[] dx;
//	dx = nullptr;
//	for (int i = 0; i < n; i++) {
//		delete[] jacobian[i];
//	}
//	delete[] jacobian;
//	jacobian = nullptr;
//	delete[] eps;
//	eps = nullptr;
//}
//
//
//void explictEuler(double* y, double** a, double* b, double f1(double**, double*, double*, int), double f2(double**, double*, double*, int)
//	, double f3(double**, double*, double*, int), double local_err, double tau_max, double T, int n) {
//	double* f = new double[n];
//	double* tau = new double[n];
//	double t = 0, min_tau = 0;
//	int it_num = 0;
//	typedef double (*funcmas_t)(double** a, double* b, double* u, int n);		//объявл. + инициал. массива функций
//	funcmas_t functions[3];
//	functions[0] = f1;
//	functions[1] = f2;
//	functions[2] = f3;										//
//
//
//	ofstream fout("Euler method.txt", ios::trunc);
//	fout << "explict Euler:\nt	u1	u2	u3:\n";
//	while (t < T) {
//		for (int i = 0; i < n; i++) {				//вычисление вектора f
//			f[i] = functions[i](a, b, y, n);
//		}											//
//		//print(f, n);
//
//		for (int i = 0; i < n; i++) {				//определение шага интегрирования
//			tau[i] = local_err / (fabs(f[i]) + (local_err / tau_max));
//		}
//		min_tau = tau[0];
//		for (int i = 0; i < n; i++) {
//			if (tau[i] < min_tau) {
//				min_tau = tau[i];
//			}
//		}						//
//
//
//		if (it_num % 10 == 0) {
//			fout << t;
//			for (int i = 0; i < n; i++) {
//				fout << "	" << y[i];
//			}
//			fout << "\n";
//		}
//
//
//
//		for (int i = 0; i < n; i++) {		//нахождение y[k+1] и t[k+1]
//			y[i] += min_tau * f[i];
//		}
//		t += min_tau;							//
//		it_num++;
//	}
//	fout.close();
//
//
//	delete[] f;
//	f = nullptr;
//	delete[] tau;
//	tau = nullptr;
//}
//
//
