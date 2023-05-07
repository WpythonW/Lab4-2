#include <iostream>
#include <cmath>
using namespace std;

struct out {
    double diff;
    double x;
    double y;
    double z;
};

// Define the function to solve
double func(double x, double y, double z) {
    if (x != 0)
        return ((2 * x + 4) * z - 2 * y) / (x * (x + 4));
    else
        return (4 * z + 2 * y) / 4;
}

// Define the derivative of the function
double deriv(double x, double y, double z) {
    return z;
}

// Define the 4th order Runge-Kutta method
out rk4(double x0, double y0, double z0, double xf, double h) {
    double x = x0, y = y0, z = z0;
    out Ret;
    //cout << "RungeKutta method" << "\n";
    while (x < xf) {
        // Calculate the k values
        double k1y = h * deriv(x, y, z);
        double k1z = h * func(x, y, z);
        double k2y = h * deriv(x + h / 2, y + k1y / 2, z + k1z / 2);
        double k2z = h * func(x + h / 2, y + k1y / 2, z + k1z / 2);
        double k3y = h * deriv(x + h / 2, y + k2y / 2, z + k2z / 2);
        double k3z = h * func(x + h / 2, y + k2y / 2, z + k2z / 2);
        double k4y = h * deriv(x + h, y + k3y, z + k3z);
        double k4z = h * func(x + h, y + k3y, z + k3z);

        // Calculate the new values of y and z
        y = y + (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
        z = z + (k1z + 2 * k2z + 2 * k3z + k4z) / 6;

        // Update the value of x
        x = x + h;

        // Output the current values of x and y
        //cout << "x: " << x << " y: " << y << " y': " << z << endl;
    }
    Ret.diff = y - z - 3;
    Ret.x = 2;
    Ret.y = y;
    Ret.z = z;
    return Ret;
}

double Fp(double x0, double y0, double z0, double xf, double h) {
    return (-rk4(x0, y0 + 2 * h, z0, xf, h).diff + 8 * rk4(x0, y0 + h, z0, xf, h).diff - 8 * rk4(x0, y0 - h, z0, xf, h).diff + rk4(x0, y0 - 2 * h, z0, xf, h).diff) / (12 * h);
    //return (rk4(x0, y0 + hs, z0, xf, h).diff - rk4(x0, y0 - hs, z0, xf, h).diff) / (2 * hs);
}

int main() {
    double x0 = 0, y0 = 5, z0 = 1.0, xf = 2.0, h = 0.01, eps = 0.01;
    out val = rk4(x0, y0, z0, xf, h);
    while (abs(val.diff) > eps) {
        y0 = y0 - val.diff / Fp(x0, y0, z0, xf, h);
        val = rk4(x0, y0, z0, xf, h);
    }
    cout << "y(0) = " << y0 << endl;
    cout << "x: " << val.x << " y: " << val.y << " y': " << val.z << endl;
    return 0;
}
