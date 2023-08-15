#include <cmath>
#include <iostream>
#include <complex>
#include <algorithm>
#include <Eigen/Dense>
#define M_PI 3.14159

std::complex<double> sine(std::complex<double> z) {;
	return std::sin(z);
}

std::complex<double> cosine(std::complex<double> z) {;
	return std::cos(z);
}

std::complex<double> j0(std::complex<double> z) {
	return std::sin(z) / (z);
}

std::complex<double> j0P(std::complex<double> z) {
	return ((z)*std::cos(z) - std::sin(z)) / pow(z,2);
}

std::complex<double> y0(std::complex<double> z) {
	return -std::cos(z) / (z);
}

std::complex<double> y0P(std::complex<double> z) {
	return ((z) * std::sin(z) - std::cos(z)) / pow(z, 2);
}

std::complex<double> linear(std::complex<double> z, std::complex<double> a, std::complex<double> b) {
	return a * z + b;
}

std::complex<double> linearP(std::complex<double> z, std::complex<double> a) {
	return a;
}


std::complex<double> quadratic(std::complex<double> z,std::complex<double> a, std::complex<double> b, std::complex<double> c) {
	return a * pow(z, 2) + b * z + c;
}

std::complex<double> quadraticP(std::complex<double> z, std::complex<double> a, std::complex<double> b) {
	return 2.*a*z + b ;
}

std::complex<double> cubic(std::complex<double> z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) {
	return a * pow(z, 3) + b * pow(z,2) + c * z + d ;
}

std::complex<double> cubicP(std::complex<double> z, std::complex<double> a, std::complex<double> b, std::complex<double> c) {
	return 3. * a * pow(z, 2) + 2. * b * z + c;
}


Eigen::MatrixXcd sine(Eigen::MatrixXcd z) {
	return z.array().sin();
}

Eigen::MatrixXcd cosine(Eigen::MatrixXcd z) {
	return z.array().cos();
}

Eigen::MatrixXcd tangent(Eigen::MatrixXcd z) {
	return z.array().tan();
}

Eigen::MatrixXcd secant(Eigen::MatrixXcd z) {
	return 2.*z.array().cos() / ((z.array().cos()* z.array().cos() - z.array().sin()* z.array().sin()) +1);
}


Eigen::MatrixXcd j0(Eigen::MatrixXcd z) {
	return z.array().sin() / z.array();
}

Eigen::MatrixXcd j0P(Eigen::MatrixXcd z) {
	return (z.array()*z.array().cos() -z.array().sin()) / (z.array()*z.array());
}

Eigen::MatrixXcd j1(Eigen::MatrixXcd z) {
	return z.array().sin() / (z.array() * z.array()) - (z.array().cos()/ z.array());
}

Eigen::MatrixXcd j1P(Eigen::MatrixXcd z) {
	return ((z.array() * z.array() - 2.) * z.array().sin() + 2 * z.array() * z.array().cos()) / (z.array() * z.array() * z.array());
}

Eigen::MatrixXcd y0(Eigen::MatrixXcd z) {
	return -z.array().cos() / z.array();
}

Eigen::MatrixXcd y0P(Eigen::MatrixXcd z) {
	return (z.array()*z.array().sin() - z.array().cos()) / (z.array() * z.array());
}



Eigen::MatrixXcd cubic(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) {
	return a * pow(z.array(), 3) + b * pow(z.array(), 2) + c * z.array() + d;
}

Eigen::MatrixXcd cubicP(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c) {
	return 3.*a * pow(z.array(), 2) + 2.*b *z.array() + c;
}

Eigen::MatrixXcd sinh(Eigen::MatrixXcd z) {
	return z.array().sinh();
}

Eigen::MatrixXcd cosh(Eigen::MatrixXcd z) {
	return z.array().cosh();
}

Eigen::MatrixXcd customDriver1(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) { //d/dz(cos(z) (z^3 - a) e^(-(cos(z) - b i))) = e^(-cos(z) + i b) ((a - z^3) sin(z) + cos(z) ((z^3 - a) sin(z) + 3 z^2))
	std::complex<double> expOffset(0., 1./M_PI);
	return z.array().cos() * (pow(z.array(), 3) - a) * exp(-(z.array().cos() - b * expOffset));
}


Eigen::MatrixXcd customDriver1P(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) { // d/dz(cos(z) (z^3 - a) e^(-(cos(z) - b i))) = e^(-cos(z) + i b) ((a - z^3) sin(z) + cos(z) ((z^3 - a) sin(z) + 3 z^2))
	std::complex<double> expOffset(0., 1./M_PI);
	return  exp(-(z.array().cos() - b * expOffset)) * (a - pow(z.array(), 3) * z.array().sin() + z.array().cos() * ((pow(z.array(), 3) - a) * z.array().sin() + 3 * pow(z.array(), 2)));
}
