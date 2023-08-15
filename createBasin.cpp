#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <complex>
#include <Config.h>
#include <driverFunctions.h>
#include <plotBasins.h>
#include <fstream>
#include <iomanip> 
#include <tuple>


void fillFinalPlane(Eigen::MatrixXcd &finalPlane, Eigen::MatrixXd& iterationMap, Eigen::MatrixXcd &dz, Eigen::MatrixXcd &plane, Config &config, int it) {
	double tol = config.getTol();
	for (int i = 0; i < dz.size(); i++) {
		if (std::abs(dz(i)) < tol && std::abs(finalPlane(i).real()) == 0 && std::abs(finalPlane(i).imag()) == 0) {
			finalPlane(i) = plane(i);
			iterationMap(i) = it;
		}
	}
}

void printIteration(int it, int iterations) {
	std::cout << "Working on iteration: " << it+1 << " of " << iterations << std::endl;
}


std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> iterateMatrixSine(const Eigen::MatrixXcd plane0, Config& config, std::complex<double> offset) {

	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();
	Eigen::MatrixXd iterationMap = Eigen::MatrixXd::Zero(plane.rows(), plane.cols());


	for (int i = 0; i < config.getIterations(); i++) {
		printIteration(i, config.getIterations());
		Eigen::MatrixXcd dz = sine(plane).array() / cosine(plane).array();
		fillFinalPlane(finalPlane, iterationMap, dz, plane, config, i);
		plane = plane.array() - dz.array();
	}
	return { finalPlane, iterationMap };
}

std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> iterateMatrixTan(const Eigen::MatrixXcd plane0, Config& config, std::complex<double> offset) {

	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();
	Eigen::MatrixXd iterationMap = Eigen::MatrixXd::Zero(plane.rows(), plane.cols());

	for (int i = 0; i < config.getIterations(); i++) {
		printIteration(i, config.getIterations());;
		Eigen::MatrixXcd dz = tangent(plane).array() / (secant(plane).array()* secant(plane).array());
		fillFinalPlane(finalPlane, iterationMap, dz, plane, config, i);
		plane = plane.array() - dz.array();
	}
	return { finalPlane, iterationMap };
}

std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> iterateMatrixSinh(const Eigen::MatrixXcd plane0, Config& config, std::complex<double> offset) {

	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();
	Eigen::MatrixXd iterationMap = Eigen::MatrixXd::Zero(plane.rows(), plane.cols());

	for (int i = 0; i < config.getIterations(); i++) {
		printIteration(i, config.getIterations());
		Eigen::MatrixXcd dz = sinh(plane).array() / cosh(plane).array();
		fillFinalPlane(finalPlane, iterationMap, dz, plane, config, i);
		plane = plane.array() - dz.array();
	}
	return { finalPlane, iterationMap };
}
std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> iterateMatrixBesselSK(const Eigen::MatrixXcd plane0, Config& config, std::complex<double> offset) {

	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();
	Eigen::MatrixXd iterationMap = Eigen::MatrixXd::Zero(plane.rows(), plane.cols());

	for (int i = 0; i < config.getIterations(); i++) {
		printIteration(i, config.getIterations());
		Eigen::MatrixXcd dz = y0(plane).array() / y0P(plane).array();
		fillFinalPlane(finalPlane, iterationMap, dz, plane, config, i);
		plane = plane.array() - dz.array();
	}
	return { finalPlane, iterationMap };
}
std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> iterateMatrixBessel(const Eigen::MatrixXcd plane0, Config& config, std::complex<double> offset ) {
	Eigen::MatrixXcd plane = plane0.array();
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();
	Eigen::MatrixXd iterationMap = Eigen::MatrixXd::Zero(plane.rows(), plane.cols());

	for (int i = 0; i < config.getIterations(); i++) {
		printIteration(i, config.getIterations());
		Eigen::MatrixXcd dz = (j0(plane).array() + offset) / j0P(plane).array();
		fillFinalPlane(finalPlane, iterationMap, dz, plane, config, i);
		plane = plane.array() - dz.array();
	}
	return { finalPlane, iterationMap };
}

std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> iterateMatrixBesselTwoTerm(const Eigen::MatrixXcd plane0, Config& config, std::complex<double> offset) {

	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();
	Eigen::MatrixXd iterationMap = Eigen::MatrixXd::Zero(plane.rows(), plane.cols());

	for (int i = 0; i < config.getIterations(); i++) {
		printIteration(i, config.getIterations());
		Eigen::MatrixXcd dz = (j0(plane).array() + j1(plane).array()) / (j0P(plane).array() + j1P(plane).array());
		fillFinalPlane(finalPlane, iterationMap, dz, plane, config, i);
		plane = plane.array() - dz.array();
	}
	return {finalPlane, iterationMap};
}


std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> iterateMatrixCubic(const Eigen::MatrixXcd plane0, Config& config, std::complex<double> offset) {

	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();
	Eigen::MatrixXd iterationMap = Eigen::MatrixXd::Zero(plane.rows(), plane.cols());

	std::complex<double> a(1, 0);
	std::complex<double> b(2, 0);
	std::complex<double> c(3, 0);
	std::complex<double> d(-1, 0);


	for (int i = 0; i < config.getIterations(); i++) {
		printIteration(i, config.getIterations());
		Eigen::MatrixXcd dz = cubic(plane, a, b, c, d).array() / cubicP(plane, a, b, c).array();
		fillFinalPlane(finalPlane, iterationMap, dz, plane, config, i);
		plane = plane.array() - dz.array();
	}
	return std::make_tuple(finalPlane, iterationMap);
}


std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> iterateMatrixCustomDriver1(const Eigen::MatrixXcd plane0, Config& config, std::complex<double> offset) {
	Eigen::MatrixXcd plane = plane0;
	Eigen::MatrixXcd finalPlane = plane0.array() - plane0.array();
	Eigen::MatrixXd iterationMap = Eigen::MatrixXd::Zero(plane.rows(), plane.cols());

	std::complex<double> a(1, 0);
	std::complex<double> b(0, 0);
	std::complex<double> c(0, 0);
	std::complex<double> d(0, -1.*M_PI);


	for (int i = 0; i < config.getIterations(); i++) {
		printIteration(i, config.getIterations());
		Eigen::MatrixXcd dz = customDriver1(plane, a, b, c, d).array() / customDriver1P(plane, a, b, c, d).array();
		fillFinalPlane(finalPlane, iterationMap, dz, plane, config, i);
		plane = plane.array() - dz.array();
	}
	return std::make_tuple(finalPlane, iterationMap);
}

void exportData(Eigen::MatrixXcd finalPlane, Eigen::MatrixXcd iterationMap, Config& config, std::complex<double> offset, int offsetIt) {
	
	std::string driverString = config.getDriver();
	std::string domainString = config.getDomainString(config);
	
	std::vector<std::vector<double>> finalPlaneVectors = {};
	double screenStepSize = (config.getScreenXmax() - config.getScreenXmin()) / config.getScreenDivs();
	double screenSize = (config.getScreenXmax() - config.getScreenXmin());

	Eigen::VectorXd finalPlaneWF(5);
	for (int i = 0; i < finalPlane.rows()* finalPlane.cols(); i++) {
		finalPlaneWF << finalPlane(i).real(), finalPlane(i).imag(), (std::floor(i / config.getScreenDivs())* screenStepSize) - screenSize/2., ((i % config.getScreenDivs()) * screenStepSize) - screenSize/2. , iterationMap(i).real();
		std::vector<double> finalPlaneWFVector(finalPlaneWF.data(), finalPlaneWF.data() + finalPlaneWF.size());
		finalPlaneVectors.push_back(finalPlaneWFVector);
	}

	std::string savePath = "C:\\Users\\Michael\\Documents\\Programming\\NewtonBasinImages\\"+ driverString +"\\variableOffset\\" + domainString + "_" + std::to_string(config.getIterations()) + "_" + std::to_string(config.getScreenDivs()) + "_" + std::to_string(offsetIt)+ "_[" + std::to_string(config.getOffsetReal()) + "," + std::to_string(config.getOffsetImag()) + "].csv";

	std::ofstream out(savePath);
	for (int i=0; i < finalPlaneVectors.size(); i++) {
		for (int j=0; j < 5; j++) {
			out << std::setprecision(8) << finalPlaneVectors[i][j] << ',';
		}
		out << '\n';
	}
}

int main() {
	int numberImages = 180;
	for (int offsetN = 0; offsetN < numberImages; offsetN++) {

		//std::complex<double> off(0.70710 * float((float(offsetN) / float(numberImages))), 0.70710 *(pow(float(offsetN)/float(numberImages),2)) );
		std::complex<double> off(std::cos((2*M_PI*float(float( (2. * offsetN) + 45.) / float(2. * numberImages)))) , std::sin((2 * M_PI * float(float((2 * offsetN) + 45.) / float(2 * numberImages)))));


		Config config(-M_PI*4., M_PI * 4., -M_PI * 4., M_PI * 4., off, "Bessel", 50, 1.e-12, 1000);
		Eigen::MatrixXcd plane = config.makeScreen(config);
		std::string domainString = config.getDomainString(config);

		if (config.getDriver() == "Cubic") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixCubic(plane, config, off);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, off, offsetN);
		}
		else if (config.getDriver() == "Sine") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixSine(plane, config, off);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, off, offsetN);
		}
		else if (config.getDriver() == "Sinh") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixSinh(plane, config, off);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, off, offsetN);
		}
		else if (config.getDriver() == "Tan") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXcd> finalOutput = iterateMatrixTan(plane, config, off);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, off, offsetN);
		}
		else if (config.getDriver() == "Bessel") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixBessel(plane, config, off);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, off, offsetN);
		}
		else if (config.getDriver() == "BesselSK") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixBesselSK(plane, config, off);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, off, offsetN);
		}
		else if (config.getDriver() == "BesselTwoTerm") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixBesselTwoTerm(plane, config, off);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, off, offsetN);
		}
		else if (config.getDriver() == "CustomDriver1") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixCustomDriver1(plane, config, off);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, off, offsetN);
		}
		else {
			std::cout << config.getDriver();
			return 1;
		}
	}
	return 0;
}