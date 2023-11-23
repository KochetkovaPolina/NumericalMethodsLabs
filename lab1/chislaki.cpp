#include <iostream>
#include <cmath>
#include <vector>


void ZeroFill(std::vector<double>& result, int& size) {
    for (int index = 0; index < size; index++)
        result.push_back(0);
}

void MultiplyMatrix(std::vector<std::vector<double> >& A, std::vector<double>& x, std::vector<double>& matrixMultiplication) {
    int matrixSize = A.size();
    ZeroFill(matrixMultiplication, matrixSize);

    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            matrixMultiplication[i] += A[i][j] * x[j];
        }
    }
}

void CalculateResVector(std::vector<std::vector<double> >& A, std::vector<double>& b, std::vector<double>& result, std::vector<double>& residualVector) {
    int matrixSize = b.size();
    std::vector<double>matrixMultiplication;
    MultiplyMatrix(A, result, matrixMultiplication);
    for (int i = 0; i < matrixSize; i++) {
        residualVector.push_back(b[i] - matrixMultiplication[i]);
    }
    return;
}

void SolveSLE(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& result) {
    int matrixSize = A.size();

    // Forward propogation 
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

    // Backward propogation 
    ZeroFill(result, matrixSize);
    for (int i = matrixSize - 1; i > -1; i--) {
        double sum = 0;
        for (int j = i; j < matrixSize; j++) {
            sum += A[i][j] * result[j];
        }
        result[i] = (b[i] - sum) / A[i][i];
    }

}


int main() {
    setlocale(LC_ALL, "Ru");
    std::vector<std::vector<double> > A = {
                            {2.30, 5.70, -0.80},
                            {3.50, -2.70, 5.30},
                            {1.70, 2.30, -1.80}
    };
    std::vector<double> b = { -6.49, 19.20, -5.09 };
    std::vector<double> result;

    std::vector<std::vector<double> > A1(A);
    std::vector<double> b1(b);
     
    std::vector<double> residualVector;

    SolveSLE(A, b, result);
    CalculateResVector(A1, b1, result, residualVector);

    std::cout <<"решение слау: ";
    for (auto iter : result)    

    std::cout << iter << " ";

    std:: cout<< std::endl;
   
    std::cout << "вектор невязки: ";
    for (auto iter : residualVector)
        std::cout << iter << " ";

    std::cout << std::endl;

    double maxvalue = 0;

    for (auto iter : residualVector) {
        if (abs(iter) > maxvalue)
        maxvalue = iter;
    }
    std::cout << "максимальное значение векторы невязки: "<< maxvalue;
    return 0;
}