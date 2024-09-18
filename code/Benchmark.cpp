#include <iostream>
#include <cmath>
#include <chrono>
using namespace std;

const double pi = 3.14159265358979323846;

void generateCartesianPoints(int n, double* x, double* y, double* z) {
    for (int i = 0; i < n; ++i) {
        x[i] = rand() % 1000;
        y[i] = rand() % 1000;
        z[i] = rand() % 1000;
    }
}

void generatePolarPoints(int n, double* r, double* theta) {
    for (int i = 0; i < n; ++i) {
        r[i] = rand() % 1000;
        theta[i] = (rand() % 360) * pi / 180.0;
    }
}

void generateSphericalPoints(int n, double* r, double* theta, double* phi) {
    for (int i = 0; i < n; ++i) {
        r[i] = rand() % 1000;
        theta[i] = (rand() % 180) * pi / 180.0;
        phi[i] = (rand() % 360) * pi / 180.0;
    }
}

double cartesianDistance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

double cartesianDistance3D(double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
}

double polarDistance(double r1, double theta1, double r2, double theta2) {
    return sqrt(pow(r1, 2) + pow(r2, 2) - 2 * r1 * r2 * cos(theta2 - theta1));
}

double sphericalDistance(double r1, double theta1, double phi1, double r2, double theta2, double phi2) {
    return sqrt(pow(r1, 2) + pow(r2, 2) - 2 * r1 * r2 * (sin(theta1) * sin(theta2) * cos(phi1 - phi2) + cos(theta1) * cos(theta2)));
}

double sphericalSurfaceDistance(double r, double theta1, double phi1, double theta2, double phi2) {
    return r * (acos(sin(phi1) * sin(phi2) + cos(phi1) * cos(phi2) * cos(theta1 - theta2)));
}

int main() {
    const int n = 100000;
    double* x = new double[n];
    double* y = new double[n];
    double* z = new double[n];

    double* r = new double[n];
    double* theta = new double[n];
    double* phi = new double[n];

    generateCartesianPoints(n, x, y, z);

    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < n - 1; ++i) {
        cartesianDistance(x[i], y[i], x[i + 1], y[i + 1]);
    }
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> cartesian2DTime = end - start;
    cout << "Calculation time for Cartesian coordinates (2D): " << cartesian2DTime.count() << " seconds" << endl;

    start = chrono::high_resolution_clock::now();
    for (int i = 0; i < n - 1; ++i) {
        cartesianDistance3D(x[i], y[i], z[i], x[i + 1], y[i + 1], z[i + 1]);
    }
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> cartesian3DTime = end - start;
    cout << "Calculation time for Cartesian coordinates (3D): " << cartesian3DTime.count() << " seconds" << endl;

    generatePolarPoints(n, r, theta);

    start = chrono::high_resolution_clock::now();
    for (int i = 0; i < n - 1; ++i) {
        polarDistance(r[i], theta[i], r[i + 1], theta[i + 1]);
    }
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> polarTime = end - start;
    cout << "Calculation time for polar coordinates:" << polarTime.count() << " seconds" << endl;

    generateSphericalPoints(n, r, theta, phi);

    start = chrono::high_resolution_clock::now();
    for (int i = 0; i < n - 1; ++i) {
        sphericalDistance(r[i], theta[i], phi[i], r[i + 1], theta[i + 1], phi[i + 1]);
    }
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> sphericalVolumeTime = end - start;
    cout << "Calculation time for spherical coordinates (by volume): " << sphericalVolumeTime.count() << " seconds" << endl;

    start = chrono::high_resolution_clock::now();
    for (int i = 0; i < n - 1; ++i) {
        sphericalSurfaceDistance(r[i], theta[i], phi[i], theta[i + 1], phi[i + 1]);
    }
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> sphericalSurfaceTime = end - start;
    cout << "Calculation time for spherical coordinates (on the surface): " << sphericalSurfaceTime.count() << " seconds" << endl;

    delete[] x;
    delete[] y;
    delete[] z;
    delete[] r;
    delete[] theta;
    delete[] phi;

    return 0;
}
