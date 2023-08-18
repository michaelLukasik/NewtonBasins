#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <complex>
#include <Config.h>
#include <driverFunctions.h>
#include <fstream>
#include <iomanip> 
#include <tuple>
#include <string>
#include "csv.h"


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
		Eigen::MatrixXcd dz = (sine(plane).array() + offset) / cosine(plane).array();
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
		Eigen::MatrixXcd dz = (tangent(plane).array() + offset) / (secant(plane).array()* secant(plane).array());
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
		Eigen::MatrixXcd dz = (sinh(plane).array() + offset) / cosh(plane).array();
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
		Eigen::MatrixXcd dz = (y0(plane).array() + offset) / y0P(plane).array();
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
		Eigen::MatrixXcd dz = (j0(plane).array() + j1(plane).array() + offset) / (j0P(plane).array() + j1P(plane).array());
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
		Eigen::MatrixXcd dz = (cubic(plane, a, b, c, d).array() + offset) / cubicP(plane, a, b, c).array();
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

int getTrackPoints(std::string fileString) {
	int count = 0;
	std::string fileLine;
	std::ifstream trackFile(fileString);
	if (trackFile.is_open()){
		while (trackFile.peek() != EOF){
			getline(trackFile, fileLine);
			count++;
		}
		trackFile.close();

		return count;
	}
	else
		std::cout << "Couldn't open the file\n";
		return 0;
}

std::complex<double> getOffsetFromTrack(std::string fileString, int rowIndex) {
	std::vector<double> xPosVector;
	std::vector<double> yPosVector;
	io::CSVReader<2> in(fileString);
	in.read_header(io::ignore_extra_column, "X", "Y");
	std::string xpos; std::string ypos;
	while (in.read_row(xpos, ypos)) {
		xPosVector.push_back(stod(xpos));
		yPosVector.push_back(stod(ypos));
	}
	std::cout << std::complex<double>(xPosVector[rowIndex], yPosVector[rowIndex]) << std::endl;
	return std::complex<double>(xPosVector[rowIndex], yPosVector[rowIndex]);
}
	

int main() {
	std::string trackFileString = "C:\\Users\\Michael\\Documents\\Programming\\NewtonBasins\\python\\tracks\\barcalona\\barcalonaTrackPoints.csv";
	std::ifstream trackFile(trackFileString);
	std::string lineInFile;
	int trackPoints = 0;
	for (int rowIndex = 0; rowIndex < getTrackPoints(trackFileString); rowIndex++) {
		std::complex<double> offset = getOffsetFromTrack(trackFileString, rowIndex);
		
		std::cout << "Offset: " << offset << std::endl;

		Config config(-M_PI*4., M_PI *4., -M_PI*4., M_PI * 4., offset, "Bessel", 50, 1.e-12, 500);
		Eigen::MatrixXcd plane = config.makeScreen(config);
		std::string domainString = config.getDomainString(config);

		if (config.getDriver() == "Cubic") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixCubic(plane, config, offset);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, offset, rowIndex);
		}
		else if (config.getDriver() == "Sine") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixSine(plane, config, offset);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, offset, rowIndex);
		}
		else if (config.getDriver() == "Sinh") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixSinh(plane, config, offset);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, offset, rowIndex);
		}
		else if (config.getDriver() == "Tan") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXcd> finalOutput = iterateMatrixTan(plane, config, offset);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, offset, rowIndex);
		}
		else if (config.getDriver() == "Bessel") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixBessel(plane, config, offset);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, offset, rowIndex);
		}
		else if (config.getDriver() == "BesselSK") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixBesselSK(plane, config, offset);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, offset, rowIndex);
		}
		else if (config.getDriver() == "BesselTwoTerm") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixBesselTwoTerm(plane, config, offset);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, offset, rowIndex);
		}
		else if (config.getDriver() == "CustomDriver1") {
			std::tuple<Eigen::MatrixXcd, Eigen::MatrixXd> finalOutput = iterateMatrixCustomDriver1(plane, config, offset);
			exportData(std::get<0>(finalOutput), std::get<1>(finalOutput), config, offset, rowIndex);
		}
		else {
			std::cout << config.getDriver();
			return 1;
		}
	}
	return 0;
}