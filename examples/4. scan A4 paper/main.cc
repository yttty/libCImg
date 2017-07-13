/*
 * File: main.cc
 * Description: detect edges of a paper and warp to a standard A4 image
 * Author: Terrill Yang
 * Date: Apr 5, 2017
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

void testA4Wraping(const string& infile) {
	//***** read image *****************//
	CImgUtil img(infile.c_str());

	array<pair<int,int>, 4> corners = img.four_corners(5, 10, 650, 200, 1, false, 25, true, false);

	CImg<float> result = img.a4_wraping(corners);
	result.display();
	string outf("result/");
	outf+=infile;
	result.save_bmp(outf.c_str());
}

int main(int argc, char const *argv[])
{
	if (argc != 2) {
		cerr << "Usage: "<< argv[0] <<" {infile path}" << endl;
		exit(-1);
	}
	else {
		testA4Wraping(argv[1]);
	}

	return 0;
}