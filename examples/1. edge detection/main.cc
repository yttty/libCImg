/*
 * File: main.cc
 * Description: Test case for class hough and class canny
 * Author: Terrill Yang
 * Date: Mar 10, 2017
 */

#include <iostream>
#include <string>
#include "canny.h"
#include "hough.h"


void test(string infile)
{
	string outfileGradient = "Gradient_Pro.bmp"; /* saving the (normalised) gradient to file? */
	string outfileNMS = "Edge_Output.bmp"; /* saving the binary canny edges to file? */
    canny cny(infile, outfileGradient, outfileNMS, 1.5f, 6.0f);

    infile = "Edge_Output.bmp";
	string outfile = "Hough_Output.bmp";
	hough hog(infile, outfile, 0.0f, 1.0f, 200.0f, 0.4f);
}

int main(int argc, char** argv)
{
	if (argc != 2) 
	{
		printf("Usage: %s {input filename}\n", argv[0]);
		return 1;
	}	

    test(argv[1]);

    return 0;
}

