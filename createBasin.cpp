#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

#include <Config.h>
#include <driverFunctions.h>
#include <plotBasins.h>
#include <fstream>
#include <iomanip> 

#define M_PI 3.14159612





//std::complex<double> iterate(std::complex<double> z0, std::complex<double> (*fcnPtr)(std::complex<double>, std::complex<double>), std::complex<double> (*drvPtr)(std::complex<double>, std::complex<double>)) {
//	return z0;
//}
Eigen::MatrixXcd iterateMatrixSine(const Eigen::MatrixXcd plane0, Config& config) {

	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();

	for (int i = 0; i < config.getIterations(); i++) {
		std::cout << "Iteration: " << i << std::endl;
		Eigen::MatrixXcd dz = sine(plane).array() / cosine(plane).array();
		for (int i = 0; i < dz.size(); i++) {
			if (std::abs(dz(i)) < config.getTol() && std::abs(finalPlane(i).real()) == 0 && std::abs(finalPlane(i).imag()) == 0) {
				finalPlane(i) = plane(i);
			}
		}
		plane = plane.array() - dz.array();
	}
	return finalPlane;
}

Eigen::MatrixXcd iterateMatrixTan(const Eigen::MatrixXcd plane0, Config& config) {

	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();

	for (int i = 0; i < config.getIterations(); i++) {
		std::cout << "Iteration: " << i << std::endl;
		Eigen::MatrixXcd dz = tangent(plane).array() / (secant(plane).array()* secant(plane).array());
		for (int i = 0; i < dz.size(); i++) {
			if (std::abs(dz(i)) < config.getTol() && std::abs(finalPlane(i).real()) == 0 && std::abs(finalPlane(i).imag()) == 0) {
				finalPlane(i) = plane(i);
			}
		}
		plane = plane.array() - dz.array();
	}
	return finalPlane;
}

Eigen::MatrixXcd iterateMatrixSinh(const Eigen::MatrixXcd plane0, Config& config) {

	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();

	for (int i = 0; i < config.getIterations(); i++) {
		std::cout << "Iteration: " << i << std::endl;
		Eigen::MatrixXcd dz = sinh(plane).array() / cosh(plane).array();
		for (int i = 0; i < dz.size(); i++) {
			if (std::abs(dz(i)) < config.getTol() && std::abs(finalPlane(i).real()) == 0 && std::abs(finalPlane(i).imag()) == 0) {
				finalPlane(i) = plane(i);
			}
		}
		plane = plane.array() - dz.array();
	}
	return finalPlane;
}
Eigen::MatrixXcd iterateMatrixBesselSK(const Eigen::MatrixXcd plane0, Config& config) {

	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();

	for (int i = 0; i < config.getIterations(); i++) {
		std::cout << "Iteration: " << i << std::endl;
		Eigen::MatrixXcd dz = y0(plane).array() / y0P(plane).array();
		for (int i = 0; i < dz.size(); i++) {
			if (std::abs(dz(i)) < config.getTol() && std::abs(finalPlane(i).real()) == 0 && std::abs(finalPlane(i).imag()) == 0) {
				finalPlane(i) = plane(i);
			}
		}
		plane = plane.array() - dz.array();
	}
	return finalPlane;
}

Eigen::MatrixXcd iterateMatrixBessel(const Eigen::MatrixXcd plane0, Config& config) {
	
	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();
	
	for (int i = 0; i < config.getIterations(); i++) {
		std::cout << "Iteration: " << i << std::endl;
		Eigen::MatrixXcd dz = j0(plane).array() / j0P(plane).array();
		for (int i = 0; i < dz.size(); i++) {
			if (std::abs(dz(i)) < config.getTol() && std::abs(finalPlane(i).real()) == 0 && std::abs(finalPlane(i).imag()) == 0){
				finalPlane(i) = plane(i);
			}
		}
		plane = plane.array() - dz.array();
	}
	return finalPlane;
}


Eigen::MatrixXcd iterateMatrixCubic(const Eigen::MatrixXcd plane0, Config& config) {

	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();
	std::complex<double> a(1, 0);
	std::complex<double> b(-0.3, 0);
	std::complex<double> c( 0.22, 0);
	std::complex<double> d(4.5, 0);
	double tol = config.getTol();

	for (int i = 0; i < config.getIterations(); i++) {
		std::cout << "Iteration: " << i+1 << " of " << config.getIterations() <<  std::endl;
		Eigen::MatrixXcd dz = cubic(plane,a,b,c,d).array() / cubicP(plane,a,b,c).array();
		for (int i = 0; i < dz.size(); i++) {
			if (std::abs(dz(i)) < tol && std::abs(finalPlane(i).real()) == 0 && std::abs(finalPlane(i).imag()) == 0) {
				finalPlane(i) = plane(i);
			}
		}
		plane = plane.array() - dz.array();
	}
	return finalPlane;
}

void exportData(Eigen::MatrixXcd finalPlane, Config& config) {
	std::vector<std::vector<double>> finalPlaneVectors = {};
	double screenStepSize = (config.getScreenXmax() - config.getScreenXmin()) / config.getScreenDivs();
	double screenSize = (config.getScreenXmax() - config.getScreenXmin());
	Eigen::VectorXd finalPlaneWF(4);
	for (int i = 0; i < finalPlane.rows()* finalPlane.cols(); i++) {
		finalPlaneWF << finalPlane(i).real(), finalPlane(i).imag(), (std::floor(i / config.getScreenDivs())* screenStepSize) - screenSize/2., ((i % config.getScreenDivs()) * screenStepSize) - screenSize/2.;
		std::vector<double> finalPlaneWFVector(finalPlaneWF.data(), finalPlaneWF.data() + finalPlaneWF.size());
		finalPlaneVectors.push_back(finalPlaneWFVector);
	}

	std::string savePath = "C:\\Users\\Michael\\Documents\\Programming\\NewtonBasinImages\\Tan\\1.csv";

	std::ofstream out(savePath);
	for (int i=0; i < finalPlaneVectors.size(); i++) {
		for (int j=0; j < 4; j++) {
			out << std::setprecision(12) << finalPlaneVectors[i][j] << ',';
		}
		out << '\n';
	}
}

int main() {
	Config config(-M_PI, M_PI , -M_PI , M_PI , 15, 1e-5, 1000);
	Eigen::MatrixXcd plane = config.makeScreen(config);
	config.coutParams();


	//Eigen::MatrixXcd finalPlane = iterateMatrixBesselSK(plane, config);
	//Eigen::MatrixXcd finalPlane = iterateMatrixBessel(plane, config);
	//Eigen::MatrixXcd finalPlane = iterateMatrixCubic(plane, config);
	//Eigen::MatrixXcd finalPlane = iterateMatrixSine(plane, config);
	//Eigen::MatrixXcd finalPlane = iterateMatrixSinh(plane, config);
	Eigen::MatrixXcd finalPlane = iterateMatrixTan(plane, config);

	
	exportData(finalPlane, config);
	//std::cout << finalPlane << std::endl;
	return 0;
}