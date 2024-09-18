#include <iostream>
#include <cmath>
using namespace std;

const double pi = 3.14159265358979323846;

void cartesianToPolar(double x, double y, double& r, double& theta) {
    r = sqrt(pow(x, 2) + pow(y, 2));
    theta = atan2(y, x); 
}

void polarToCartesian(double r, double theta, double& x, double& y) {
    x = r * cos(theta);  
    y = r * sin(theta); 
}

int main() {
    double cartesian_x[] = { 5.0, 10.0, 2.0 };
    double cartesian_y[] = { 4.0, 5.0, 3.0 };
    int num_points = 3;

    for (int i = 0; i < num_points; i++) {
        double r, theta;
        cartesianToPolar(cartesian_x[i], cartesian_y[i], r, theta);
        cout << "Cartesian coordinates (" << cartesian_x[i] << ", " << cartesian_y[i] << ") in the polar system: (r = " << r << ", theta = " << theta << ")" << endl;
    }

    cout << endl;

    double polar_r[] = { 6.0, 11.0, 3.6 };
    double polar_theta[] = { pi / 5, pi / 6, pi / 3 };

    for (int i = 0; i < num_points; i++) {
        double x, y;
        polarToCartesian(polar_r[i], polar_theta[i], x, y);
        cout << "Polar coordinates (r = " << polar_r[i] << ", theta = " << polar_theta[i] << ") in the Cartesian system: (" << x << ", " << y << ")" << endl;
    }

    return 0;
}
