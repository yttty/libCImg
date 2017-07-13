/*
 * File: main.cc
 * Description: Define test code for histeq
 * Author: Terrill Yang
 * Date: Mar 15, 2017
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <string.h>
#include "CImg.h"
#include "libcimg.h"
using namespace cimg_library;
using namespace std;

void testHistEq(const string& infile) {
	//***** read image *****************//
	CImgUtil img(infile.c_str()); 

	CImg<float> res1 = img.histeq_gray();
	CImg<float> res2 = img.histeq_rgb();
	CImg<float> res3 = img.histeq_hsi();
	(img.get_img_rgb(), res1).display();
	(img.get_img_rgb(), res2).display();
	(img.get_img_rgb(), res3).display();
	res1.save_bmp("histeq_gray.bmp");
	res2.save_bmp("histeq_rgb.bmp");
	res3.save_bmp("histeq_hsi.bmp");
}

int main(int argc, char const *argv[])
{
	if (argc != 2) {
		cerr << "Usage: "<< argv[0] <<" {infile path}" << endl;
		exit(-1);
	} 
	else {
		// testHistEq
		testHistEq(argv[1]);
	}

	return 0;
}