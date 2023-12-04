#include <iostream>
#include <iomanip>

using namespace std;

double* gauss(double** matrix, int n) {
	double* X = new double[n];
	for (int i = 0; i < n; i++) {
		double max = fabs(matrix[i][i]);
		int save = i;
		for (int j = i + 1; j < n; j++) {
			if (matrix[j][i] > max) {
				max = fabs(matrix[i][j]);
				save = j;
			}
		}
		if (fabs(max) < 1e-6) {
			cout << "Решения нет!" << endl;
			return 0;
		}
		swap(matrix[i], matrix[save]);
		double div = matrix[i][i];
		for (int j = 0; j < n + 1; j++) {
			matrix[i][j] /= div;
		}
		for (int j = i + 1; j < n; j++) {
			for (int k = n; k >= i; k--) {
				matrix[j][k] -= matrix[j][i] * matrix[i][k];
			}
		}
	}
	for (int i = n - 1; i >= 0; i--) {
		X[i] = matrix[i][n];
		for (int j = 0; j < i; j++) {
			matrix[j][n] -= matrix[j][i] * matrix[i][n];
		}
	}
	return X;
}

void LS(double* x, double* y, int N, int m) {
	double* POWERX = new double[2 * m];
	for (int i = 0; i < 2 * m; i++) {
		POWERX[i] = 0;
		for (int j = 0; j < N; j++) {
			POWERX[i] += pow(x[j], i + 1);
		}
	}
	double** SUMX = new double* [m + 1];
	for (int i = 0; i < m + 1; i++) {
		SUMX[i] = new double[m + 2];
	}
	for (int i = 0; i < m + 1; i++) {
		for (int j = 0; j < m + 1; j++) {
			if (i + j) {
				SUMX[i][j] = POWERX[i + j - 1];
			}
			else {
				SUMX[i][j] = N;
			}
		}
	}

	/*for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < m; ++j)
			cout << SUMX[i][j] << "   ";
		cout << endl;
	}*/
	cout << endl;
	double* PRAW = new double[m + 1];
	for (int i = 0; i < m + 1; i++) {
		PRAW[i] = 0;
		for (int j = 0; j < N; j++) {
			PRAW[i] += y[j] * pow(x[j], i);
		}
	}
	for (int i = 0; i < m + 1; i++) {
		SUMX[i][m + 1] = PRAW[i];
	}
	double* a = gauss(SUMX, m + 1);
	double S2 = 0;
	for (int i = 0; i < N; i++) {
		double sum = y[i];
		for (int j = 0; j < m + 1; j++) {
			sum -= a[j] * pow(x[i], j);
		}
		S2 += pow(sum, 2);
	}
	S2 /= double(N - m - 1);
	double sigma = sqrt(S2);
	cout << "Коэфициенты a: " << endl;
	for (int i = 0; i < m + 1; i++) {
		cout << a[i] << " ";
	}
	cout << "\nCреднеквадратическое отклонение: " << sigma << endl;
	cout << "Полином: y = ";
	for (int i = m; i > 0; i--) {
		if (i == 1) {
			cout << a[i] << "*x" << " + " << a[i - 1] << endl;
		}
		else cout << a[i] << "*x^" << i << " + ";
	}
	delete[] POWERX;
	for (int i = 0; i < m + 1; i++) {
		delete[] SUMX[i];
	}
	delete[] SUMX;
	delete[] PRAW;
	delete[] a;
}

int main() {
	system("chcp 1251");
	system("cls");
	const int N = 5;
	const int m = 1;
	double* x = new double [N] {5.76, 12.25, 25, 47.4721, 100}; //вариант 10, х возведенны в квадрат из-за замены
	double* y = new double [N] {0.0141, 0.0281, 0.0562, 0.1125, 0.2250};
	cout << setw(2) << "x" << "\t" << "y" << endl;
	for (int i = 0; i < N; i++) {
		cout << setw(2) << x[i] << "\t" << y[i] << endl;
	}
	LS(x, y, N, m);
	delete[] x;
	delete[] y;
}