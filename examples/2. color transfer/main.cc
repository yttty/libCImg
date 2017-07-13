/*
 * File: main.cc
 * Description: Define test code for color transfer
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

void testColorTransfer(const string& infile, const string& targetfile) {
	//***** read image *****************//
	CImgUtil img(infile.c_str()); 
	CImg<float> target(targetfile.c_str()); 

	CImg<float> res = img.color_transfer(target);
	(img.get_img_rgb(), target, res).display();
	res.save_bmp("color_transfer.bmp");
}

int main(int argc, char const *argv[])
{
	if (argc != 2) {
		cerr << "Usage: "<< argv[0] <<" {target image path} {result path}" << endl;
		exit(-1);
	} 
	else {
		// testColortransfer
		testColorTransfer(argv[1], argv[2]);
	}

	return 0;
}