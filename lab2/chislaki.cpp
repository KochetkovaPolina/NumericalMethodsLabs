#pragma once
#include <iostream>
#include <cmath>
#include <vector>


void Zeros(std::vector<double>& result, int& size) {
    for (int index = 0; index < size; index++)
        result.push_back(0);
}

void MatrixMultiplication(std::vector<std::vector<double> >& A, std::vector<double>& x, std::vector<double>& matrixMultiplication) {
    int matrixSize = A.size();
    Zeros(matrixMultiplication, matrixSize);

    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            matrixMultiplication[i] += A[i][j] * x[j];
        }
    }

    return;
}

void ResidualVectorCalculation(std::vector<std::vector<double> >& A, std::vector<double>& b, std::vector<double>& result, std::vector<double>& residualVector) {
    int matrixSize = b.size();
    std::vector<double>matrixMultiplication;
    MatrixMultiplication(A, result, matrixMultiplication);
    for (int i = 0; i < matrixSize; i++) {
        residualVector.push_back(b[i] - matrixMultiplication[i]);
    }
    return;
}

void GaussElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& result) {
    int matrixSize = A.size();
    for (int i = 0; i < matrixSize; i++) {
        int indexMax = i;
        for (int j = i + 1; j < matrixSize; j++) {
            if (abs(A[j][i]) > abs(A[indexMax][i])) {
                indexMax = j;
            }
        }

        std::swap(A[i], A[indexMax]);
        std::swap(b[i], b[indexMax]);

        for (int j = i + 1; j < matrixSize; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < matrixSize; k++) {
                A[j][k] -= factor * A[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    Zeros(result, matrixSize);
    for (int i = matrixSize - 1; i > -1; i--) {
        double sum = 0;
        for (int j = i; j < matrixSize; j++) {
            sum += A[i][j] * result[j];
        }
        result[i] = (b[i] - sum) / A[i][i];
    }
    return;
}


std::vector<double> solve(std::vector<std::vector<double>> A, std::vector<double> b) {
    setlocale(LC_ALL, "Ru");
    
    std::vector<double> result;

    std::vector<std::vector<double> > A1(A);
    std::vector<double> b1(b);

    std::vector<double> residualVector;

    GaussElimination(A, b, result);
    ResidualVectorCalculation(A1, b1, result, residualVector);

    return result;
    
}