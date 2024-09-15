﻿#include <iostream>
#include <cmath>
using namespace std;

const double pi = 3.14159265358979323846;

double degreesToRadians(double degrees) {
    return degrees * pi / 180.0;
}

double distance2D(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

double distance3D(double x1, double y1, double z1, double x2, double y2, double z2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
}

double polarDistance(double r1, double theta1, double r2, double theta2) {
    return sqrt(pow(r1, 2) + pow(r2, 2) - 2 * r1 * r2 * cos(degreesToRadians(theta2 - theta1)));
}

void sphericalToCartesian(double r, double theta, double phi, double& x, double& y, double& z) {
    x = r * sin(degreesToRadians(phi)) * cos(degreesToRadians(theta));
    y = r * sin(degreesToRadians(phi)) * sin(degreesToRadians(theta));
    z = r * cos(degreesToRadians(phi));
}

double sphericalDistance(double r1, double theta1, double phi1, double r2, double theta2, double phi2) {
    double x1, y1, z1;
    double x2, y2, z2;
    sphericalToCartesian(r1, theta1, phi1, x1, y1, z1);
    sphericalToCartesian(r2, theta2, phi2, x2, y2, z2);
    return distance3D(x1, y1, z1, x2, y2, z2);
}

double sphericalSurfaceDistance(double r, double theta1, double phi1, double theta2, double phi2) {
    double deltaSigma = acos(sin(degreesToRadians(phi1)) * sin(degreesToRadians(phi2)) + cos(degreesToRadians(phi1)) * cos(degreesToRadians(phi2)) * cos(degreesToRadians(theta2 - theta1)));
    return r * deltaSigma;
}

int main() {
    double x1 = 1, y1 = 2, x2 = 4, y2 = 6;
    cout << "2D Distance: " << distance2D(x1, y1, x2, y2) << endl;

    double z1 = 3, z2 = 7;
    cout << "3D Distance: " << distance3D(x1, y1, z1, x2, y2, z2) << endl;

    double r1 = 5, theta1 = 30, r2 = 10, theta2 = 60;
    cout << "Polar Distance: " << polarDistance(r1, theta1, r2, theta2) << endl;

    double phi1 = 45, phi2 = 60;
    cout << "Spherical Distance: " << sphericalDistance(r1, theta1, phi1, r2, theta2, phi2) << endl;

    double sphereRadius = 10;
    cout << "Spherical Surface Distance: " << sphericalSurfaceDistance(sphereRadius, theta1, phi1, theta2, phi2) << endl;

    return 0;
}
