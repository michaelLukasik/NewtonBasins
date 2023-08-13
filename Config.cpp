#include <iostream>
#include <cmath>
#include <Eigen/Dense>


#include "Config.h"

Eigen::MatrixXcd Config::makeScreen(Config& config) {
	double screenXSpacing = (config.getScreenXmax() - config.getScreenXmin()) / config.getScreenDivs();
	double screenYSpacing = (config.getScreenYmax() - config.getScreenYmin()) / config.getScreenDivs(); // unused
	double screenLength = config.getScreenXmax() - config.getScreenXmin();

	Eigen::MatrixXcd complexPlane(config.getScreenDivs(), config.getScreenDivs());

	for (int realPixel = 0; realPixel < config.getScreenDivs(); realPixel++) {
		for (int imPixel = 0; imPixel < config.getScreenDivs(); imPixel++) {
			
			std::complex<double> val( (realPixel*screenXSpacing) - screenLength/2., (imPixel * screenYSpacing) - screenLength / 2.);
			
			complexPlane(realPixel, imPixel) = val;
		}
	}

	return complexPlane;
}