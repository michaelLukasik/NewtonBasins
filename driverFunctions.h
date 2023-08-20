#pragma once
#include <complex>

std::complex<double> sine(std::complex<double> z);
std::complex<double> cosine(std::complex<double> z);

std::complex<double> j0(std::complex<double> z);
std::complex<double> j0P(std::complex<double> z);

std::complex<double> y0(std::complex<double> z);
std::complex<double> y0P(std::complex<double> z);

std::complex<double> linear(std::complex<double> z, std::complex<double> a, std::complex<double> b);
std::complex<double> linearP(std::complex<double> z, std::complex<double> a);

std::complex<double> quadratic(std::complex<double> z, std::complex<double> a, std::complex<double> b, std::complex<double> c);
std::complex<double> quadraticP(std::complex<double> z, std::complex<double> a, std::complex<double> b);

std::complex<double> cubic(std::complex<double> z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d);
std::complex<double> cubicP(std::complex<double> z, std::complex<double> a, std::complex<double> b, std::complex<double> c);



Eigen::MatrixXcd sine(Eigen::MatrixXcd z);
Eigen::MatrixXcd cosine(Eigen::MatrixXcd z);
Eigen::MatrixXcd tangent(Eigen::MatrixXcd z);
Eigen::MatrixXcd secant(Eigen::MatrixXcd z);

Eigen::MatrixXcd sinh(Eigen::MatrixXcd z);
Eigen::MatrixXcd cosh(Eigen::MatrixXcd z);

Eigen::MatrixXcd j0(Eigen::MatrixXcd z);
Eigen::MatrixXcd j0P(Eigen::MatrixXcd z);

Eigen::MatrixXcd j1(Eigen::MatrixXcd z);
Eigen::MatrixXcd j1P(Eigen::MatrixXcd z);

Eigen::MatrixXcd y0(Eigen::MatrixXcd z);
Eigen::MatrixXcd y0P(Eigen::MatrixXcd z);

Eigen::MatrixXcd cubic(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d);
Eigen::MatrixXcd cubicP(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c);

Eigen::MatrixXcd zeta(Eigen::MatrixXcd z);
Eigen::MatrixXcd zetaP(Eigen::MatrixXcd z);


Eigen::MatrixXcd customDriver1(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d);
Eigen::MatrixXcd customDriver1P(Eigen::MatrixXcd z, std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d);

void zetaPTest();