#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double P(double x) {
    if (x != 0)
        //return -x;
        return -(2.0 * x + 4.0) / (x * (x + 4.0));
    else
        return 0.5; //взял lim x->0
}

double Q(double x) {
    //return -1;
    return 2 / (x * (x + 4));
}

double f(double x) {
    return 0;
}

void prVect(std::vector<double> X) {
    for (double x : X) {
        std::cout << x << ' ';
    }
    std::cout << std::endl;
}

void prMatr(std::vector<vector<double>> M) {
    for (std::vector<double> row : M) {
        for (double x : row) {
            std::cout << x << " ";
        }
        cout << std::endl;
    }
}

vector<double> GenerateB(int N, std::vector<double> X, double h) {
    std::vector<double> B(N);
    //B[0] = h + (h * h)/2 * -P(X[0]);
    //B[0] = -0.98;
    B[0] = h;
    for (int i = 1; i < N - 1; i++) {
        B[i] = h * h * f(X[i]);
    }
    B[N - 1] = 18 * h;

    //B[N - 1] = 0;
    return B;
}


vector<vector<double>> Generate3DiagMatrix(int N, std::vector<double> X, double h) {
    std::vector<vector<double>> Matr(N, std::vector<double>(3));

    //Matr[0][0] = 0;
    //Matr[0][1] = -1 + h * h * Q(X[0]) - 1/2 * P(X[0]);
    //Matr[0][2] = 1 + 1/2 * P(X[0]);

    Matr[0][0] = 0.0;
    Matr[0][1] = -1.0;
    Matr[0][2] = 1.0;
    // 
    //Matr[0][0] = 0;
    //Matr[0][1] = -2.04;
    //Matr[0][2] = 1.02;

    for (int row = 1; row < N - 1; row++) {
        Matr[row][0] = 1.0 - P(X[row]) * h / 2;
        Matr[row][1] = -2.0 + h * h * Q(X[row]);
        Matr[row][2] = 1.0 + P(X[row]) * h / 2.0;
        //cout << X[row] << endl;
    }

    Matr[N - 1][0] = -2.0; 
    Matr[N - 1][1] = h + 2.0; 
    Matr[N - 1][2] = 0.0;

    //Matr[N - 1][0] = 1;
    //Matr[N - 1][1] = -1.4;
    //Matr[N - 1][2] = 0;

    return Matr;
}


std::vector<double> tridiagonalSolve(const std::vector<double>& b, const std::vector<std::vector<double>>& Matr) {
    int n = b.size();
    std::vector<double> alpha(n - 1);
    std::vector<double> beta(n);
    std::vector<double> y(n);

    // Проверка условия диагонального преобладания
    for (int i = 0; i < n; i++) {
        double sum = abs(Matr[i][0]) + abs(Matr[i][1]);
        if (abs(Matr[i][2]) > sum) {
            cout << Matr[i][0] << " " << Matr[i][1] << " " << Matr[i][2] << endl;
            std::cerr << "Matrix doesn't meet conditions of diagonal supersion" << std::endl;
            return y;
        }
    }

    alpha[0] = -Matr[0][2] / Matr[0][1];
    beta[0] = b[0] / Matr[0][1];

    for (int i = 1; i < n - 1; i++) {
        double denom = Matr[i][1] + Matr[i][0] * alpha[i - 1];
        alpha[i] = -Matr[i][2] / denom;
        beta[i] = (b[i] - Matr[i][0] * beta[i - 1]) / denom;
    }

    beta[n - 1] = (b[n - 1] - Matr[n - 1][0] * beta[n - 2]) / (Matr[n - 1][1] + Matr[n - 1][0] * alpha[n - 2]);

    y[n - 1] = beta[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        y[i] = alpha[i] * y[i + 1] + beta[i];
    }

    return y;
}

int main() {
    double start = 0.0;
    double end = 2.0;
    double h = 0.0001;
    int N = (end - start) / h;

    std::vector<double> X;
    for (double i = start + h; i < end + h; i += h) {
        X.push_back(i);
    }

    std::vector<double> B = GenerateB(N, X, h);
    std::vector<vector<double>> Matr = Generate3DiagMatrix(N, X, h);

    //prVect(X);
    //prMatr(Matr);
    //prVect(B);

    std::vector<double> Y = tridiagonalSolve(B, Matr);

    cout << "Y values: " << endl;
    prVect(Y);

    return 0;
}
