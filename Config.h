#pragma once

#include <iostream>
#include <complex>
#include <string>

#define M_PI 3.14159



class Config
{
public: 
	Config(double screenXmin, double screenXmax, double screenYmin, double screenYmax, std::complex<double> offset, std::string driver ,int iterations = 1, double tol = 0.001, int screenDivs = 100)
		: m_screenXmin{ screenXmin }, m_screenXmax{ screenXmax }, m_screenYmin{ screenYmin }, m_screenYmax{ screenYmax }, m_offset{ offset },
		m_driver{ driver }, m_iterations{ iterations }, m_tol{ tol }, m_screenDivs{ screenDivs }
	{
	}
	void coutParams() {
		std::cout << m_screenXmin << " : screen Xmin" << std::endl;
		std::cout << m_screenXmax << " : screen Xmax" << std::endl;
		std::cout << m_screenYmin << " : screen Ymin" << std::endl;
		std::cout << m_screenYmax << " : screen Ymax" << std::endl;
		std::cout << m_iterations << " : iterations" << std::endl;
		std::cout << m_tol << " : tolerance" << std::endl;
		std::cout << m_driver << " : Driver Function" << std::endl;
	}
	
	typedef void (*functionCall)(int args);
	typedef void (*derivativeCall)(int args);

	double getScreenXmin() {return  m_screenXmin;}
	double getScreenXmax() { return  m_screenXmax; }
	double getScreenYmin() { return  m_screenYmin; }
	double getScreenYmax() { return  m_screenYmax; }

	
	int getIterations() { return m_iterations; }
	int getScreenDivs() { return  m_screenDivs; }

	std::string getDriver() { return m_driver; }
	
	double getTol() { return m_tol; }
	
	double getOffsetReal() { return m_offset.real(); }
	double getOffsetImag() { return m_offset.imag(); }


	Eigen::MatrixXcd makeScreen(Config config);

	std::string getDomainString(Config &config);



private: 

	double m_screenXmin;
	double m_screenXmax;
	double m_screenYmin;
	double m_screenYmax;

	std::string m_driver;

	int m_screenDivs;
	int m_iterations;
	
	double m_tol;
	
	std::complex<double> m_offset =(1.,1.);

};

