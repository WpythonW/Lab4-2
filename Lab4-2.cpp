#include <iostream>
#include <cmath>
using namespace std;

int main() {
    double h = 0.1;
    int n = 21;
    double* x = new double[n];
    double* y = new double[n];
    double** A = new double* [n];
    for (int i = 0; i < n; i++) {
        A[i] = new double[n];
        x[i] = i * h;
        y[i] = 0;
    }

    // Начальное условие y'(0) = 1
    y[1] = 1;

    // Граничное условие y(2) - y'(2) = 3
    double alpha = 1 / h;
    double beta = -1 / h - 3;

    // Заполнение матрицы коэффициентов и правой части уравнения
    double* b = new double[n];

    for (int i = 1; i < n - 1; i++) {
        A[i][i - 1] = -1 / pow(h, 2) - (x[i] + 4) / 2 / h;
        A[i][i] = 2 / pow(h, 2) + x[i] * (x[i] + 4);
        A[i][i + 1] = -1 / pow(h, 2) + (x[i] - 4) / 2 / h;
        b[i] = 0;
    }

    // Граничные условия
    A[0][0] = 1;
    b[0] = y[0];

    A[n - 1][n - 2] = alpha;
    A[n - 1][n - 1] = beta;
    b[n - 1] = -alpha * y[n - 2];

    // Решение системы уравнений методом прогонки
    for (int i = 1; i < n; i++) {
        double m = A[i][i - 1] / A[i - 1][i - 1];
        A[i][i] -= m * A[i - 1][i];
        b[i] -= m * b[i - 1];
    }

    y[n - 1] = b[n - 1] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        y[i] = (b[i] - A[i][i + 1] * y[i + 1]) / A[i][i];
    }

    // Вывод результата
    for (int i = 0; i < n; i++) {
        cout << "y(" << x[i] << ") = " << y[i] << endl;
    }

    // Освобождение памяти
    delete[] x;
    delete[] y;
    for (int i = 0; i <
        n; i++) {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;

    return 0;
}
