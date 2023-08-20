#include <cmath>
#include <iostream>
#include <complex>
#include <algorithm>
#include <Eigen/Dense>
#include <numbers>

#include <vector>


#define M_PI 3.14159

const double LOWER_THRESHOLD = 1.0e-6;
const double UPPER_BOUND = 1.0e+4;
const int maxTermsInZetaSummation = 50;


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

Eigen::MatrixXcd zeta(Eigen::MatrixXcd z) { // Eigen only supports exponents >1 and shifts >0 https://eigen.tuxfamily.org/dox/unsupported/namespaceEigen.html

	Eigen::MatrixXcd sum = Eigen::MatrixXcd::Zero(z.rows(), z.cols());
	Eigen::MatrixXcd initialArray = Eigen::MatrixXcd::Zero(z.rows(), z.cols());
	Eigen::MatrixXcd sumFlattened = sum.reshaped<Eigen::RowMajor>().transpose();
	Eigen::MatrixXcd zFlattened = z.reshaped<Eigen::RowMajor>().transpose();

	Eigen::MatrixXcd terms = Eigen::MatrixXcd::Zero(maxTermsInZetaSummation, zFlattened.cols());


	std::complex<double> half(0.5, 0.0);
	std::complex<double> one(1.0, 0.0);
	std::complex<double> two(2.0, 0.0);
	std::complex<double> rev(-1.0, 0.0);


	
	terms.row(0) = half / (one - pow(two, (one - zFlattened.array())));
	for (int n = 1; n < maxTermsInZetaSummation; n++){
		std::complex<double> nCplx(n, 0.0); //complex index
		for (int k = 0; k < n; k++){
			std::complex<double> kCplx(k, 0.0); //complex index
			terms.row(k) *= half * (nCplx / (nCplx - kCplx));
			sumFlattened += terms.row(k);
		}
		terms.row(n).array() = (rev * terms.row(n - 1).array() * pow((nCplx / (nCplx + one)), zFlattened.array()) / nCplx);
		sumFlattened.array() += terms.row(n).array();
	}
	return sumFlattened.reshaped(z.rows(), z.cols());
}


Eigen::MatrixXcd zetaP(Eigen::MatrixXcd z){
	std::complex<double> one(1.0, 0.0);
	std::complex<double> im(0.0, 1.0);

	Eigen::MatrixXcd sum = Eigen::MatrixXcd::Zero(z.rows(), z.cols());
	for (int i = 1; i < maxTermsInZetaSummation; i++) {
		sum.array() += (log(double(i)) / pow(i, z.array()));
	}
	return -1. * sum;
}


void zetaPTest()
{
	Eigen::MatrixXcd z(2,2);
	std::complex<double> im(0.0, 1.0);
	z(0, 0) = 0.73 + 9.11 * im;
	z(0, 1) = 0.1 + 1.7 * im;
	z(1, 0) = 10.9 - 7. * im;
	z(1, 1) = 0.01 - 0.2 * im;
	std::complex<double> testingZ(0.73, 9.11);
	std::complex<double> one(1.0, 0.0);
	std::complex<double> sumTestingZ(0.,0.);

	std::complex<double> testZ(1, 0);
	Eigen::MatrixXcd testSum(2,2);



	for (int k = 1; k < 100; k++) {
		if (k > 98) { std::cout <<  z << std::endl; }
		testSum.array() += (  log(double(k)) / pow(k, z.array())  );
	}

	for (int k = 1; k < 100; k++) {

		sumTestingZ += std::log(double(k)) / std::pow(k, testingZ)  ; ;
		//sumTestingZ += testingZ;
	}
}


