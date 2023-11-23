#include "chislaki.cpp"
#include <iostream>
#include <vector>
#include <cmath>


// Variant 3
// x1^2 * x2^2 - 3x1^2 - 6x2^3 + 8 = 0
// x1^4 - 9x2 + 2 = 0
// starts (-1.5, 1.5) and (-1, 1)



double func1(double x1, double x2) {
    return x1 * x1 * x2 * x2 - 3 * x1 * x1 - 6 * x2 * x2 * x2 + 8;
}

double func1_dx1(double x1, double x2) {
    return 2 * x1 * x2 * x2 - 6 * x1;
}

double func1_dx2(double x1, double x2) {
    return 2 * x1 * x1 * x2 - 18 * x2 * x2;
}

double func2(double x1, double x2) {
    return x1 * x1 * x1 * x1 - 9 * x2 + 2;
}

double func2_dx1(double x1, double x2) {
    return 4 * x1 * x1 * x1;
}

double func2_dx2(double x1, double x2) {
    return -9;
}


void hand_made(double x1, double x2, double E1, double E2, int max_iter) {
    double delta1 = std::max(fabs(func1(x1, x2)),
                             fabs(func2(x1, x2)));
    double delta2 = 1;

    int iter = 0;

    while ((delta1 > E1 || delta2 > E2) && iter < max_iter) {
        iter++;

        std::cout << iter << ":\t" << "x1: " << x1 << "\t" << "x2: " << x2 << std::endl;

        std::vector<double> minus_F{-func1(x1, x2), -func2(x1, x2)};
        std::vector<std::vector<double>> J{{func1_dx1(x1, x2), func1_dx2(x1, x2)},
                                           {func2_dx1(x1, x2), func2_dx2(x1, x2)}};

        std::vector<double> solution = solve(J, minus_F);
        x1 += solution[0];
        x2 += solution[1];

        delta1 = fabs(minus_F[0]);
        for (int i = 1; i < minus_F.size(); i++) {
            if (delta1 < fabs(minus_F[i]))
                delta1 = fabs(minus_F[i]);
        }

        double v1 = fabs(x1) < 1 ? solution[0] : solution[0] / x1;
        double v2 = fabs(x2) < 1 ? solution[1] : solution[1] / x2;;
        delta2 = std::max(v1, v2);
    }
}

std::vector<std::vector<double>> get_J(double x1, double x2, double k) {
    return std::vector<std::vector<double>> {{(func1(x1 + x1 * k, x2) - func1(x1, x2)) / k / x1, (func1(x1, x2 + x2 * k) - func1(x1, x2)) / k / x2},
                                             {(func2(x1 + x1 * k, x2) - func2(x1, x2)) / k / x1, (func2(x1, x2 + x2 * k) - func2(x1, x2)) / k / x2}};

}

void ebanina3000(double x1, double x2, double E0, double E1, int max_iter, double k) {
    double delta0 = std::max(func1(x1, x2),
                             func2(x1, x2));
    double delta1 = 1;

    int iter = 0;

    while ((delta0 > E0 || delta1 > E1) && iter < max_iter) {
        iter++;
        std::cout << iter << ":\t" << "x1: " << x1 << "\t" << "x2: " << x2 << std::endl;

        std::vector<double> F{-func1(x1, x2), -func2(x1, x2)};
        std::vector<std::vector<double>> J = get_J(x1, x2, k);

        std::vector<double> solution = solve(J, F);
        x1 += solution[0];
        x2 += solution[1];

        delta0 = fabs(F[0]);
        for (int i = 1; i < F.size(); i++) {
            if (delta0 < fabs(F[i]))
                delta0 = fabs(F[i]);
        }

        double v1 = fabs(x1) < 1 ? solution[0] : solution[0] / x1;
        double v2 = fabs(x2) < 1 ? solution[1] : solution[1] / x2;;
        delta1 = std::max(v1, v2);
    }
}



int main() {

    double x1 = -1.5;
    double x2 = 1.5;

    double E0 = 1e-5;
    double E1 = 1e-5;
    int max_iter = 100;

    hand_made(x1, x2, E0, E1, max_iter);
    std::cout << '\n';
    ebanina3000(x1, x2, E0, E1, max_iter, 0.05);


    return 0;
}
