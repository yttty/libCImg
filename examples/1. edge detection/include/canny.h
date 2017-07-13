/*
 * File: canny.h
 * Description: Define class canny
 * Author: Terrill Yang
 * Date: Mar 10, 2017
 */

#ifndef _CANNY_
#define _CANNY_

#include <cmath>
#include <string>
#include <iostream>
#include "CImg.h" /* C++ image library */
using namespace cimg_library;
using namespace std;


class canny {
private:
	/* canny parameters */
	float	sigma;
	float	threshold;

private:
    void gauss_filter( CImg<float> & filter, float sigma, int deriv);

public:

    canny(string infile, string	outfileGradient, string	outfileNMS, \
    	float sigma, float threshold); //Constructor
    void CannyDiscrete( CImg<float> in, float sigma, float threshold,
         CImg<float> &outSmooth, CImg<float> &outGradient,
         CImg<float> &outOrientation, CImg<float> &outThreshold,
         CImg<float> &outNMS );
};

#endif
