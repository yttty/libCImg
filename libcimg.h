/*
 * File: main.cc
 * Description: Define test code for histeq and color transfer
 * Author: Terrill Yang
 * Date: Apr 5, 2017
 */

#ifndef LIBCIMG_H
#define LIBCIMG_H

#include "CImg.h"
#include <iostream>
#include <string>
#include <array>
#include <algorithm>
#include <string.h>

using namespace std;
using namespace cimg_library;

#define PI 3.14159265358979323846

// An encapsulation of a bunch od utilities of CImg library
class CImgUtil
{
private:
	CImg<float> img_rgb_;
	int width_, height_;

	// 2D Dot(x, y) [value]
	struct Dot {
	    int x, y, value;
	    Dot(int _x, int _y, int _value)
	        :x(_x), y(_y), value(_value){}
	};

	// 2D Straight Line(x, y) [y = m*x + b]
	struct Line {
	    float m, b;
	    Line(float _m, float _b)
	        :m(_m), b(_b){}
	};

private:
	// Do histogram equalization on one dimension(depth)
	// Usage Example: histeq_one(in_img, out_img, dimension)
	void histeq_one(const CImg<float>&, CImg<float>&, int);
	// Do histogram equalization on in tensity level(for HSI color space only)
	// Usage Example: histeq_intensity(in_img, out_img)
	// Note: in_img must be in HSI color space
	void histeq_intensity(const CImg<float>&, CImg<float>&);

	// Calculate mean and standard derivative of image
	// Usage Example: get_mean_std(img, mean, std)
	// Note: image must have three channels
	// Note: in-arg mean and std will be overwritten
	void get_mean_std(const CImg<float>&, array<float, 3>&, array<float, 3>&);

	// Gaussian Lowpass filter on gradient
	// Usage Example:
	//   CImg<float> result = non_maximum_suppression(response, threshold);
	CImg<float> non_maximum_suppression(const CImg<float>&, \
										const array<float, 3>&);

	// Gaussian Lowpass filter on gradient
	// Usage Example:
	//   CImg<float> response = gaussian_lowpass(Ixx, Ixy, Iyy, alpha, threshold);
	// Note: in-arg threshold will be overwritten
	CImg<float> gaussian_lowpass(const CImg<float>& , \
								 const CImg<float>& , \
								 const CImg<float>& , \
								 const float, array<float, 3>& );

	// Calculate the Euclidean distance
	// Usage Example: dist = distance(a, b);
	float distance(float x, float y);

	// Polar coordinate intersection at x
	const int CrossX(int theta, int distance, int x);

	// Polar coordinate intersection at y
	const int CrossY(int theta, int distance, int y);


public:
	// Constructor
	// Usage Example: CImgUtil img("filename.bmp")
	CImgUtil(const char *);


	// Do histogram equalization on gray image (color image will be
	// transfered to gray image)
	// Usage Example:
	//   CImgUtil img("filename.bmp");
	//   CImg<float> img2 = img.histeq_gray();
	CImg<float> histeq_gray();
	// Do histogram equalization on RGB channels respectively
	// Usage Example:
	//   CImgUtil img("filename.bmp");
	//   CImg<float> img2 = img.histeq_rgb();
	CImg<float> histeq_rgb();
	// Do histogram equalization on Intensity channel of HSI image
	// Usage Example:
	//   CImgUtil img("filename.bmp");
	//   CImg<float> img2 = img.histeq_hsi();
	CImg<float> histeq_hsi();


	// Do histogram equalization on Intensity channel of HSI image
	// Usage Example:
	//   CImg<float> target("targetname.bmp");
	//   CImgUtil img("filename.bmp");
	//   CImg<float> res = img.color_transfer(target);
	CImg<float> color_transfer(const CImg<float>&);

	// Do Sobel edge detection
	// Usage Example:
	// CImg<float> edge = img.sobel_edge_detect();
	CImg<float> sobel_edge_detect(unsigned int threshold = 85, \
								  unsigned int dilate = 0, \
								  bool display = false);

	// Do Harris Corner Detection
	// Usage Example:
	// CImg<float> corner = img.harris_corner_detect(0.04);
	CImg<float> harris_corner_detect(float alpha, bool heuristic = false, \
									 bool display = false);

	// Detect four corners of a paper
	// Usage Example:
	// array<pair<int,int>, 4> corners = img.four_corners(3, 20, 800, 200, 1);
	array<pair<int,int>, 4> four_corners(unsigned int blur = 5, \
										 unsigned int grad_limit = 15, \
										 unsigned int threshold = 800, \
										 unsigned int diff = 50, \
										 unsigned int slope_flag = 1, \
										 bool heuristic = true, \
										 unsigned int thres_gap = 25, \
										 bool display = false, \
										 bool save = false);

	// Detect four corners of a paper (overload)
	// Usage Example:
	// array<pair<int,int>, 4> corners = img.four_corners(3, 20, 800, 200, 1);
	array<pair<int,int>, 4> four_corners(bool display = false, \
										 bool save = false);

	// Correct A4 perspective tranmsform using the detected
	//   four corners of a paper
	// Usage Example:
	// array<pair<int,int>, 4> corners = img.four_corners(3, 20, 650, 50, 1, false, 25, false, false);
	// CImg<float> result = img.a4_wraping(corners);
	CImg<float> a4_wraping(array<pair<int,int>, 4>);

	// Get a new instance of image in RGB space
	// Usage Example:
	//   CImg<float> new = img.get_img_rgb();
	CImg<float> get_img_rgb();
	// Get a new instance of image in grayscale
	// Usage Example:
	//   CImg<float> new = img.get_img_gray();
	CImg<float> get_img_gray();
	// Get width of image
	// Usage Example:
	//   int w = img.get_width();
	int get_width();
	// Get height of image
	// Usage Example:
	//   int h = img.get_height();
	int get_height();
};

CImg<float> pointDetect(CImg<float> src);

#endif // LIBCIMG_H