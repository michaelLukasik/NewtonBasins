#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include "Config.h"



Eigen::MatrixXcd Config::makeScreen(Config config) {
	double screenXSpacing = (config.getScreenXmax() - config.getScreenXmin()) / config.getScreenDivs();
	double screenYSpacing = (config.getScreenYmax() - config.getScreenYmin()) / config.getScreenDivs(); // unused
	double screenLength = config.getScreenXmax() - config.getScreenXmin();

	Eigen::MatrixXcd complexPlane(config.getScreenDivs(), config.getScreenDivs());

	for (int realPixel = 0; realPixel < config.getScreenDivs(); realPixel++) {
		for (int imPixel = 0; imPixel < config.getScreenDivs(); imPixel++) {
			std::complex<double> val( (realPixel*screenXSpacing) - screenLength/2., (imPixel * screenYSpacing) - screenLength / 2.);
			complexPlane(realPixel, imPixel) = val + std::complex<double> (config.getOffsetReal(), config.getOffsetImag());
		}
	}

	return complexPlane;
}

std::string Config::getDomainString(Config& config) {
	double xmin = config.getScreenXmin();
	double xmax = config.getScreenXmax();
	double ymin = config.getScreenYmin();
	double ymax = config.getScreenYmax();
	double offsetRe = config.getOffsetReal();
	double offsetIm = config.getOffsetImag();



	std::string domainString = "[";
	std::vector<double> domain{ xmin,xmax,ymin,ymax };
	for (auto i : domain) {
		domainString.append(std::to_string(i) + "_");
	}
	domainString.append("]");
	std::cout << "Domain String: " << domainString << std::endl;
	std::cout << "Offset:  (" << std::to_string(offsetRe) + ", " + std::to_string(offsetIm) +")" << std::endl;
	return domainString;

}

