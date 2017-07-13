/*
 * File: hough.h
 * Description: Define class hough
 * Author: Terrill Yang
 * Date: Mar 10, 2017
 */

#ifndef _HOUGH_
#define _HOUGH_

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include "CImg.h" /* C++ image library */
using namespace cimg_library;
using namespace std;


/** a vector of pixel coordinates.
 * Usage:
 * unsigned i;
 * int x, y;
 * TVectorOfPairs nonmax;
 * nonmax.push_back (make_pair(x,y)); // adding new pixel coordinates:
 * x = nonmax[i].first; // get x-coordinate of i-th pixel
 * y = nonmax[i].second; // get y-coordinate of i-th pixel
 */
typedef std::vector<std::pair<int, int> > TVectorOfPairs;


class hough {
private:
	/* pre-scaling/rotation operations */
    float   rotate;
    float   zoom;
    /* post processing stuf */
    float   in_thresh;       /* absolute threshold for gradient magnitude */
    float   out_thresh;         /* relative to global hough space maximum */

private:
    void non_maximum_suppression( CImg<float> & img, TVectorOfPairs & nonmax,
			      float thresh, int halfwidth );
    CImg<float> overlay(CImg<float> imga, CImg<float> imgb);

public:

	hough(string, string, float, float, float, float);
	void houghTransform( const CImg<float> & img,
	    CImg<float> & HoughSpace, CImg<float> & result,
	    float in_thresh, float out_thresh );
};

#endif