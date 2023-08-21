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
const int MAXNUM = 1e2;
const int maxTermsInZetaSummation = 100;


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

Eigen::MatrixXcd tangentP(Eigen::MatrixXcd z) {
	return secant(z);
}

Eigen::MatrixXcd tangentPP(Eigen::MatrixXcd z) {
	return (2 * z.array().sin()) / (pow(z.array().cos(), 3));
}

Eigen::MatrixXcd j0(Eigen::MatrixXcd z) {
	return z.array().sin() / z.array();
}

Eigen::MatrixXcd j0P(Eigen::MatrixXcd z) {
	return (z.array()*z.array().cos() -z.array().sin()) / (z.array()*z.array());
}

Eigen::MatrixXcd j0PP(Eigen::MatrixXcd z) {
	return ((z.array() * z.array() - 2.) * z.array().sin() + 2 * z.array() * z.array().cos()) / pow(z.array(), 3);
}

Eigen::MatrixXcd j1(Eigen::MatrixXcd z) {
	return z.array().sin() / (z.array() * z.array()) - (z.array().cos()/ z.array());
}

Eigen::MatrixXcd j1P(Eigen::MatrixXcd z) {
	return ((z.array() * z.array() - 2.) * z.array().sin() + 2 * z.array() * z.array().cos()) / (z.array() * z.array() * z.array());
}

Eigen::MatrixXcd j1PP(Eigen::MatrixXcd z) {
	return ((z.array() * (pow(z.array(), 2) - 6.0) * z.array().cos()) - 3.0 * (pow(z.array(), 2) - 2.0) * z.array().sin()) / pow(z.array(), 4);
}

Eigen::MatrixXcd y0(Eigen::MatrixXcd z) {
	return -z.array().cos() / z.array();
}

Eigen::MatrixXcd y0P(Eigen::MatrixXcd z) {
	return (z.array()*z.array().sin() - z.array().cos()) / (z.array() * z.array());
}

Eigen::MatrixXcd y0PP(Eigen::MatrixXcd z) {
	return ((z.array() * z.array() - 2.) * z.array().cos() - 2. * z.array() * z.array().sin()) / pow(z.array(), 3);
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

Eigen::MatrixXcd coshSinc(Eigen::MatrixXcd z) {
	return z.array().cosh() / z.array();
}

Eigen::MatrixXcd coshSincP(Eigen::MatrixXcd z) {
	return (z.array().sinh() * z.array() - z.array().cosh()) / pow(z.array(), 2);
}
Eigen::MatrixXcd coshSincPP(Eigen::MatrixXcd z) {
	return ((pow(z.array(), 2) + 2.) * z.array().cosh() - 2. * z.array() * z.array().sinh()) / pow(z.array(), 3);
}



Eigen::MatrixXcd customDriver1(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) { //d/dz(cos(z) (z^3 - a) e^(-(cos(z) - b i))) = e^(-cos(z) + i b) ((a - z^3) sin(z) + cos(z) ((z^3 - a) sin(z) + 3 z^2))
	std::complex<double> expOffset(0., 1./M_PI);
	return z.array().cos() * (pow(z.array(), 3) - a) * exp(-(z.array().cos() - b * expOffset));
}


Eigen::MatrixXcd customDriver1P(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) { // d/dz(cos(z) (z^3 - a) e^(-(cos(z) - b i))) = e^(-cos(z) + i b) ((a - z^3) sin(z) + cos(z) ((z^3 - a) sin(z) + 3 z^2))
	std::complex<double> expOffset(0., 1./M_PI);
	return  exp(-(z.array().cos() - b * expOffset)) * (a - pow(z.array(), 3) * z.array().sin() + z.array().cos() * ((pow(z.array(), 3) - a) * z.array().sin() + 3 * pow(z.array(), 2)));
}

Eigen::MatrixXcd customDriver1PP(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) { // d/dz(cos(z) (z^3 - a) e^(-(cos(z) - b i))) = e^(-cos(z) + i b) ((a - z^3) sin(z) + cos(z) ((z^3 - a) sin(z) + 3 z^2))
	std::complex<double> expOffset(0., 1. / M_PI);
	Eigen::MatrixXcd exponent = exp((-1.0 * z.array().cos() + b * expOffset));
	Eigen::MatrixXcd zacube = (pow(z.array(), 3) - a);
	Eigen::MatrixXcd term1 = zacube.array() * z.array().cos() * -1.0 * exponent.array();
	Eigen::MatrixXcd term2 = z.array().cos() * (zacube.array() * (z.array().cos() * exponent.array() + pow(z.array().sin(), 2) * exponent.array()));
	Eigen::MatrixXcd term3 = z.array().cos() * 6.0 * pow(z.array(), 2) * z.array().sin() * exponent.array();
	Eigen::MatrixXcd term4 = z.array().cos() * 6.0 * z.array() * exponent.array();
	Eigen::MatrixXcd term5 = 2.0 * z.array().sin() * (zacube.array() * z.array().sin() * exponent.array() + 3.0 * pow(z.array(), 2) * exponent.array());
	return term1.array() + term2.array() + term3.array() + term4.array() + term5.array();
}

Eigen::MatrixXcd customDriver2(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) { //f(z) == (abs(z) * sin^3(z)) / e^(z/i)
	std::complex<double> expOffset(0., 1. / M_PI);
	Eigen::MatrixXcd eiz = exp(expOffset * z.array());
	return (abs(z.array()) * pow(z.array().sin(), 3)) / eiz.array();
}

Eigen::MatrixXcd customDriver2P(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) { 
	std::complex<double> eye(0., 1.);
	Eigen::MatrixXcd eiz = exp(-1.*eye * z.array());
	return -eye * eiz.array() * abs(z.array()) * pow(z.array().sin(), 3) + ((eiz.array() * z.array() * pow(z.array().sin(), 3)) / abs(z.array())) + (3.0 * eiz.array() * abs(z.array()) * pow(z.array().sin(), 2) * z.array().cos());
}

Eigen::MatrixXcd customDriver2PP(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) { //f(z) == (abs(z) * sin^3(z)) / e^(z/i)
	std::complex<double> eye(0., 1.);
	Eigen::MatrixXcd eiz = exp(-1. * eye * z.array());
	Eigen::MatrixXcd term1 = pow(z.array().sin(), 3) * ((-2.0 * eye * eiz.array() * z.array()) / abs(z.array()) - (eiz.array() * abs(z.array())));
	Eigen::MatrixXcd term2 = eiz.array() * abs(z.array()) * (6.0 * z.array().sin() * pow(z.array().cos(), 2) - 3.0 * pow(z.array().sin(), 3));
	Eigen::MatrixXcd term3 = 6.0 * ((eiz.array() * z.array() / abs(z.array())) - eye * eiz.array() * abs(z.array())) * z.array().cos() * pow(z.array().sin(), 2);
	return term1 + term2 + term3;
}


Eigen::MatrixXcd customDriver3(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) { 
	return 3. * z.array().imag().cos() + 7. * z.array().real().sin();
}

Eigen::MatrixXcd customDriver3P(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) {
	return -3.* z.array().imag().sin() + 7. * z.array().real().cos();
}

Eigen::MatrixXcd customDriver3PP(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) { 
	return  -3. * z.array().imag().cos() + -7. * z.array().real().sin();
}

std::complex<double> computeZeta(std::complex<double> s) {
	std::complex<double> a_arr[MAXNUM + 1];
	std::complex<double> half(0.5, 0.0);
	std::complex<double> one(1.0, 0.0);
	std::complex<double> two(2.0, 0.0);
	std::complex<double> rev(-1.0, 0.0);
	std::complex<double> sum(0.0, 0.0);
	std::complex<double> prev(1.0e+20, 0.0);

	a_arr[0] = half / (one - std::pow(two, (one - s))); //initialize with a_0 = 0.5 / (1 - 2^(1-s))
	sum += a_arr[0];

	for (int n = 1; n <= MAXNUM; n++)
	{
		std::complex<double> nCplx(n, 0.0); //complex index

		for (int k = 0; k < n; k++)
		{
			std::complex<double> kCplx(k, 0.0); //complex index

			a_arr[k] *= half * (nCplx / (nCplx - kCplx));
			sum += a_arr[k];
		}

		a_arr[n] = (rev * a_arr[n - 1] * std::pow((nCplx / (nCplx + one)), s) / nCplx);
		sum += a_arr[n];


		if (std::abs(prev - sum) < LOWER_THRESHOLD)//If the difference is less than or equal to the threshold value, it is considered to be convergent and the calculation is terminated.
			break;

		if (std::abs(sum) > UPPER_BOUND)//doesn't work for large values, so it gets terminated when it exceeds UPPER_BOUND
			break;

		prev = sum;
	}
	return sum;
}
Eigen::MatrixXcd zeta(Eigen::MatrixXcd z) { ;

	Eigen::MatrixXcd sum = Eigen::MatrixXcd::Zero(z.rows(), z.cols());
	Eigen::MatrixXcd initialArray = Eigen::MatrixXcd::Zero(z.rows(), z.cols());
	Eigen::MatrixXcd sumFlattened = sum.reshaped<Eigen::RowMajor>().transpose();
	Eigen::MatrixXcd zFlattened = z.reshaped<Eigen::RowMajor>().transpose();


	Eigen::MatrixXcd terms = Eigen::MatrixXcd::Zero(maxTermsInZetaSummation, zFlattened.cols());


	std::complex<double> half(0.5, 0.0);
	std::complex<double> one(1.0, 0.0);
	std::complex<double> two(2.0, 0.0);
	std::complex<double> rev(-1.0, 0.0);

	for (int i = 0; i < z.rows(); i++) {
		for (int j = 0; j < z.rows(); j++) {
			sum(i, j) = computeZeta(z(i, j));
		}
	}
	//terms.row(0) = half / (one - pow(two, (one - zFlattened.array())));
	//sumFlattened.array() += terms.row(0).array();
	//for (int n = 1; n < maxTermsInZetaSummation; n++){
	//	std::complex<double> nCplx(n, 0.0); //complex index
	//	for (int k = 0; k < n; k++){
	//		std::complex<double> kCplx(k, 0.0); //complex index
	//		terms.row(k) *= half * (nCplx / (nCplx - kCplx));
	//		sumFlattened += terms.row(k);
	//	}
	//	terms.row(n).array() = (rev * terms.row(n - 1.).array() * pow((nCplx / (nCplx + one)), zFlattened.array()) / nCplx);
	//	sumFlattened.array() += terms.row(n).array();
	//} 
	return sum;
}


Eigen::MatrixXcd zetaP(Eigen::MatrixXcd z){
	std::complex<double> one(1.0, 0.0);
	std::complex<double> im(0.0, 1.0);
	std::complex<double> ztemp(-10.0, -10.0);
	std::complex<double> sumtemp(0.,0.);

	Eigen::MatrixXcd sum = Eigen::MatrixXcd::Zero(z.rows(), z.cols());
	for (int i = 1; i < maxTermsInZetaSummation; i++) {
		sumtemp += -1.*one*(log(double(i)) / pow(i, ztemp));
		sum.array() += -1. * one * (log(double(i)) / pow(i, z.array()));
	}
	return -1. * sum;
}


void zetaPTest()
{
	Eigen::MatrixXcd z(2,2);
	std::complex<double> im(0.0, 1.0);
	z(0, 0) = -10. + -10. * im;
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
		testSum.array() += -1.* one *(  log(double(k)) / pow(k, z.array())  );
	}

	for (int k = 1; k < 100; k++) {

		sumTestingZ += std::log(double(k)) / std::pow(k, testingZ)  ; ;
	}
}