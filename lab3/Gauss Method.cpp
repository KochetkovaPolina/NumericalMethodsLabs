#include<iostream>
#include<iomanip>
using namespace std;

void gauss(double** matrixA, int N, double* dy) {
	double tmp;
	for (int i = 0; i < N; i++)
	{
		double max = fabs(matrixA[i][i]);
		int save = i;
		for (int j = i + 1; j < N; j++) {
			if (fabs(matrixA[j][i]) > max) {
				max = fabs(matrixA[i][j]);
				save = j;
			}
		}
		std::swap(matrixA[i], matrixA[save]);

		tmp = matrixA[i][i];
		for (int j = N; j >= i; j--)
			matrixA[i][j] /= tmp;
		for (int j = i + 1; j < N; j++)
		{
			tmp = matrixA[j][i];
			for (int k = N; k >= i; k--)
				matrixA[j][k] -= tmp * matrixA[i][k];
		}
	}
	dy[N - 1] = matrixA[N - 1][N];
	for (int i = N - 2; i >= 0; i--)
	{
		dy[i] = matrixA[i][N];
		for (int j = i + 1; j < N; j++) {
			dy[i] -= matrixA[i][j] * dy[j];
		}
	}
}
