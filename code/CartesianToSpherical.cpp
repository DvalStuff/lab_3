#include <iostream>
#include <cmath>
using namespace std;

const double pi = 3.14159265358979323846;

void cartesianToSpherical(double x, double y, double z, double& r, double& theta, double& phi) {
    r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    theta = acos(z / r);                 
    phi = atan2(y, x);                    
}

void sphericalToCartesian(double r, double theta, double phi, double& x, double& y, double& z) {
    x = r * sin(phi) * cos(theta);         
    y = r * sin(phi) * sin(theta);         
    z = r * cos(phi);                 
}

int main() {
    double cartesian_x[] = { 5.0, 1.0, -2.0 };
    double cartesian_y[] = { 4.0, 5.0, -3.0 };
    double cartesian_z[] = { 6.0, 2.0, 4.0 };
    int num_points = 3;

    for (int i = 0; i < num_points; i++) {
        double r, theta, phi;
        cartesianToSpherical(cartesian_x[i], cartesian_y[i], cartesian_z[i], r, theta, phi);
        cout << "Cartesian coordinates (" << cartesian_x[i] << ", " << cartesian_y[i] << ", " << cartesian_z[i] << ") in the spherical system: (r = " << r << ", theta = " << theta << ", phi = " << phi << ")" << endl;
    }

    cout << endl;

    double spherical_r[] = { 9.0, 5.0, 5.0 };
    double spherical_theta[] = { pi / 4, pi / 3, pi / 4 };
    double spherical_phi[] = { pi / 5, pi / 3, -((2*pi) / 3) };

    for (int i = 0; i < num_points; i++) {
        double x, y, z;
        sphericalToCartesian(spherical_r[i], spherical_phi[i], spherical_theta[i], x, y, z);
        cout << "Spherical coordinates (r = " << spherical_r[i] << ", theta = " << spherical_theta[i] << ", phi = " << spherical_phi[i] << ") in the Cartesian system: (" << x << ", " << y << ", " << z << ")" << endl;
    }
    return 0;
}
