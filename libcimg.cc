/*
 * File: libcimg.cc
 * Description: Define Class CImgUtil
 * Note: Refer to libcimg.h for documents
 * Author: Terrill Yang
 * Date: Apr 5, 2017
 */

#include "libcimg.h"

CImgUtil::CImgUtil(const char * infile)
{
	const CImg<float> in_img(infile);
	img_rgb_ = in_img;

	width_ = img_rgb_.width();
	height_ = img_rgb_.height();

	if ( width_ == 0 || height_ == 0 ) {
		cerr << "Error when loading input image." << endl;
		exit(-1);
	}
}

CImg<float> CImgUtil::histeq_gray()
{
	CImg<float> img_t = img_rgb_.get_RGBtoGray();
	histeq_one(img_rgb_.get_RGBtoGray(), img_t, 0);
	return img_t;
}

CImg<float> CImgUtil::histeq_rgb()
{
	CImg<float> img_t = img_rgb_;
	histeq_one(img_rgb_, img_t, 0);
	histeq_one(img_rgb_, img_t, 1);
	histeq_one(img_rgb_, img_t, 2);
	return img_t;
}

CImg<float> CImgUtil::histeq_hsi()
{
	CImg<float> img_t = img_rgb_.get_RGBtoHSI();
	histeq_intensity(img_rgb_.get_RGBtoHSI(), img_t);
	return img_t.get_HSItoRGB();
}

void CImgUtil::get_mean_std(const CImg<float>& img, \
							array<float, 3>& mean, \
							array<float, 3>& std)
{
	float factor = 1.0 / (img.width() * img.height());
	// calculate mean
	for (int dim = 0; dim < 3; ++dim) {
		mean[dim] = 0.0;
		cimg_forXY(img, x, y)
		{
			mean[dim] += factor * img(x, y, dim);
		}
	}
	// calculate std
	for (int dim = 0; dim < 3; ++dim) {
		std[dim] = 0.0;
		cimg_forXY(img, x, y)
		{
			std[dim] += factor * (img(x, y, dim) - mean[dim]) \
						* (img(x, y, dim) - mean[dim]);
		}
		std[dim] = sqrt(std[dim]);
	}
}

CImg<float> CImgUtil::color_transfer(const CImg<float>& target)
{
	CImg<float> img_t = img_rgb_.get_RGBtoLab();
	CImg<float> target_t = target.get_RGBtoLab();
	array<float, 3> img_t_mean, target_t_mean, img_t_std, target_t_std;
	get_mean_std(img_t, img_t_mean, img_t_std);
	get_mean_std(target_t, target_t_mean, target_t_std);

	for (int dim = 0; dim < 3; ++dim) {
		cimg_forXY(img_t, x, y)
		{
			img_t(x, y, dim) = (target_t_std[dim]/img_t_std[dim]) * \
							   (img_t(x, y, dim) - img_t_mean[dim]) + \
							   target_t_mean[dim];
		}
	}

	return img_t.get_LabtoRGB();
}

CImg<float> CImgUtil::sobel_edge_detect(unsigned int threshold, \
										unsigned int dilate, \
										bool display) {
    int sobelX[3][3] = { { -1,0,1 },{ -2,0,2 },{ -1,0,1 } };
    int sobelY[3][3] = { { 1,2,1 },{ 0,0,0 },{ -1,-2,-1 } };

    CImg<float> img = img_rgb_.get_RGBtoGray();
    CImg<float> img_t(width_, height_, 1, 3);
    img_t.fill(0);

    for (int i = 1; i < width_ - 1; i++) {
        for (int j = 1; j < height_ - 1; j++) {
            float sumX = 0.0;
            float sumY = 0.0;
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    sumX += img(i - 1 + row, j - 1 + col, 0) * \
                    		sobelX[row][col];
                    sumY += img(i - 1 + row, j - 1 + col, 0) * \
                    		sobelY[row][col];
                }
            }

            img_t(i, j, 0) = sqrt(sumX*sumX+sumY*sumY) > threshold ? 255 : 0;
            img_t(i, j, 1) = img_t(i, j, 0);
            img_t(i, j, 2) = img_t(i, j, 0);
        }
    }
    if (dilate > 0) img_t.dilate(3);
    if (display) img_t.display();
    return img_t;
}

CImg<float> CImgUtil::harris_corner_detect(float alpha, bool heuristic, \
										   bool display) {
	const CImg<float>& src = img_rgb_;

	// convolution mask
    int horizontal[3][3] = { { 5,0, -5 },{ 8,0,-8 }, { 5,0,-5 } };
    int vertical[3][3] = { { 5,8, 5 },{ 0,0,0 },{ -5,-8,-5 } };

    // calculate gradient in x & y axis
    CImg<float> Ix(width_, height_, 1, 3);
    CImg<float> Iy(width_, height_, 1, 3);

    for (int i = 1; i < width_; i++) {
        for (int j = 1; j < height_ - 1; j++) {
            for (int k = 0; k < 3; k++) {
                float sumX = 0.0;
                float sumY = 0.0;
                for (int row = 0; row < 3; row++) {
                    for (int col = 0; col < 3; col++) {
                        sumX += src(row + i - 1, col + j - 1, k) * \
                        		horizontal[row][col];
                        sumY += src(row + i - 1, col + j - 1, k) * \
                        		vertical[row][col];
                    }
                }
                Ix(i, j, k) = sumX;
                Iy(i, j, k) = sumY;
            }
        }
    }

    // product of Ix and Iy
	CImg<float> Ixx(width_, height_, 1, 3);
    CImg<float> Iyy(width_, height_, 1, 3);
    CImg<float> Ixy(width_, height_, 1, 3);

    for (int i = 0; i < width_; i++) {
        for (int j = 0; j < height_; j++) {
            for (int k = 0; k < 3; k++) {
                Ixx(i, j, k) = Ix(i, j, k)*Ix(i, j, k);
                Iyy(i, j, k) = Iy(i, j, k)*Iy(i, j, k);
                Ixy(i, j, k) = Ix(i, j, k)*Iy(i, j, k);
            }
        }
    }

    // gaussian lowpass filter
    array<float, 3> threshold;
    if (heuristic)
    	alpha = min(width_, height_) / 100;
    else
    	alpha = 0.04;
    CImg<float> response = gaussian_lowpass(Ixx, Ixy, Iyy, alpha, threshold);

    // thresholding response
    for (int i = 0; i < width_; i++) {
        for (int j = 0; j < height_; j++) {
            for (int k = 0; k < 3; k++) {
                if (response(i, j, k) < threshold[k]) {
                    response(i, j, k) = 0;
                }
            }
        }
    }

    CImg<float> result = non_maximum_suppression(response, threshold);

    // display corner
    if (display) {
    	CImg<float> disp = src;
    	unsigned char color[] = { 255,255,0 };
		for (int i = 0; i < disp.width(); i++) {
		    for (int j = 0; j < disp.height(); j++) {
		        if (result(i, j, 0) == 1) {
		            disp.draw_circle(i, j, 5, color);
		        }
		    }
		}
		disp.display();
	}

    return result;
}

array<pair<int,int>, 4> CImgUtil::four_corners(unsigned int blur, \
											   unsigned int grad_limit, \
											   unsigned int threshold, \
											   unsigned int diff, \
											   unsigned int slope_flag, \
											   bool heuristic, \
											   unsigned int thres_gap, \
											   bool display, \
											   bool save) {
	CImg<float> img_gray = img_rgb_.get_RGBtoGray();

    if (heuristic) {
        img_gray.blur(std::sqrt(width_*height_)*0.0012);
    }
    else if (blur > 0) {
        img_gray.blur(blur);
    }

    //if (display) {
    //    img_gray.display();
    //}

    CImg<float> gradnum(width_, height_, 1, 1, 0);
    int maxDistance = distance(width_, height_);
    CImg<float> hough(360, maxDistance, 1, 1, 0);
    // 3*3 neighbour matrix I
    CImg_3x3(I, float);
    // calculate grad and hough transform
    if (heuristic) grad_limit = 0;
    cimg_for3x3(img_gray, x, y, 0, 0, I, float) {
        // const float ix = (2*(Inc-Ipc)+(Inp-Ipp)+(Inn-Ipn))/4;
        // const float iy = (2*(Icn-Icp)+(Inn-Inp)+(Ipn-Ipp))/4;
        const float ix = Inc-Ipc;
        const float iy = Icn-Icp;
        float grad = std::sqrt(ix*ix + iy*iy);
        gradnum(x, y) = grad;
        if (heuristic) grad_limit = (grad_limit>grad)? grad_limit:grad;
    }

    if (heuristic) grad_limit *= 0.25;
    cimg_for3x3(img_gray, x, y, 0, 0, I, float) {
        if (gradnum(x, y) > grad_limit) {
            cimg_forX(hough, angle) {
                float rangle = (float)angle*PI / 180.0;
                int polar = (int)(x*cos(rangle) + y*sin(rangle));
                if (polar >= 0 && polar < hough.height()) {
                    hough(angle, polar) += 1;
                }
            }
        }
    }

    // Find peaks
    std::vector<Dot*> peaks;
    if (heuristic) {
        threshold = (int)0.2*distance(width_, height_);
        diff = threshold/16;
    }
    do {
    	// Peak
	    cimg_forXY(hough, angle, polar) {
	        if (hough(angle, polar) > threshold) { // is peak
	            bool flag = false;
	            const int ymin = 0;
	            const int ymax = height_ - 1;
	            const int x0 = CrossY(angle, polar, ymin);
	            const int x1 = CrossY(angle, polar, ymax);

	            const int xmin = 0;
	            const int xmax = width_ - 1;
	            const int y0 = CrossX(angle, polar, xmin);
	            const int y1 = CrossX(angle, polar, xmax);

	            if ((x0 >= 0 && x0 <= xmax) || (x1 >= 0 && x1 <= xmax) || (y0 >= 0 && y0 <= ymax) || (y1 >= 0 && y1 <= ymax)) {
	                // find the nearest peak
	                for (int i = 0; i < peaks.size(); ++i) {
	                    if (distance(peaks[i]->x - angle, peaks[i]->y - polar) < diff) {
	                        flag = true;
	                        // if larger than nearest peak, then update peak
	                        if (peaks[i]->value < hough(angle, polar)) {
	                            peaks[i] = new Dot(angle, polar, hough(angle, polar));
	                        }
	                    }
	                }
	                // no nearest peak found, then add new peak
	                if (flag == false) {
	                    peaks.push_back(new Dot(angle, polar, hough(angle, polar)));
	                }
	            }
	        }
	    }
	    if (heuristic) {
	    	if (peaks.size() > 4)
	    		threshold += thres_gap;
	    	else if (peaks.size() < 4)
	    		threshold -= thres_gap;
	    }
    } while (heuristic && peaks.size() != 4);


    cout << "lines: " << endl;
    // Transform polar coordinates to rectangular coordinates
    vector<Line*> lines;
    for (int i = 0; i < peaks.size(); ++i) {
        float angle = (float)peaks[i]->x*PI / 180.0;
        float m = -cos(angle) / sin(angle);
        float b = (float)peaks[i]->y / sin(angle);
        bool flag = true;
        for (auto& l : lines) {
            if (abs(l->m - m) < 1e-3 && abs(l->b - b) < 1e-3) {
                flag = false;
                break;
            }
        }
        if (flag) {
            lines.push_back(new Line(m, b));
            cout << "y = (" << m << ") x + (" << b << ")" << endl;
        }

    }

    cout << endl << "intersections: " << endl;

    // Calculate line intersections
    vector<Dot*> intersections;
    for (int i = 0; i < lines.size(); ++i) {
        for (int j = i + 1; j < lines.size(); ++j) {
            float m0 = lines[i]->m;
            float m1 = lines[j]->m;
            float b0 = lines[i]->b;
            float b1 = lines[j]->b;
            float x = (b1 - b0) / (m0 - m1);
            float y = (m0*b1 - m1*b0) / (m0 - m1);
            if (x >= 0 && x < width_ && y >= 0 && y < height_) {
                intersections.push_back(new Dot(x, y, 0));
                cout << "(" << x << ", " << y << ")" << endl;
            }
        }
    }
    cout << endl;

    if (display) {
    	CImg<float> result(img_rgb_);
    	// Draw lines
	    for (int i = 0; i < lines.size(); ++i) {
	        const int ymin = 0;
	        const int ymax = result.height() - 1;
	        const int x0 = (float)(ymin - lines[i]->b) / lines[i]->m;
	        const int x1 = (float)(ymax - lines[i]->b) / lines[i]->m;

	        const int xmin = 0;
	        const int xmax = result.width() - 1;
	        const int y0 = xmin*lines[i]->m + lines[i]->b;
	        const int y1 = xmax*lines[i]->m + lines[i]->b;

	        const float color[] = { 255, 255, 0 };

	        if (std::abs(lines[i]->m) > slope_flag) {
	            result.draw_line(x0, ymin, x1, ymax, color);
	        }
	        else {
	            result.draw_line(xmin, y0, xmax, y1, color);
	        }
	    }

	    // Draw intersections
	    for (int i = 0; i < intersections.size(); ++i) {
	        const float color[] = { 255, 255, 0 };
	        result.draw_circle(intersections[i]->x, intersections[i]->y, 10, color);
	    }

	    // Display
	    gradnum.display();
	    result.display();
	    if (save) {
	    	gradnum.save_bmp("gradnum.bmp");
	    	result.save_bmp("result.bmp");
		}
	    // hough.display();
    }

    // sort corners
	array<pair<int,int>, 4> corners;
    float dist;
    int idx;

    dist = maxDistance;
    for (int i = 0; i < intersections.size(); ++i) {
        if (distance(intersections[i]->x - width_, intersections[i]->y) < dist) {
            dist = distance(intersections[i]->x - width_, intersections[i]->y);
            corners[1] = make_pair(intersections[i]->x, intersections[i]->y);
            idx = i;
        }
    }
    intersections.erase(intersections.begin()+idx);

    dist = maxDistance;
    for (int i = 0; i < intersections.size(); ++i) {
        if (distance(intersections[i]->x, intersections[i]->y) < dist) {
            dist = distance(intersections[i]->x, intersections[i]->y);
            corners[0] = make_pair(intersections[i]->x, intersections[i]->y);
            idx = i;
        }
    }
    intersections.erase(intersections.begin()+idx);

    dist = maxDistance;
    for (int i = 0; i < intersections.size(); ++i) {
        if (distance(intersections[i]->x, intersections[i]->y - height_) < dist) {
            dist = distance(intersections[i]->x, intersections[i]->y - height_);
            corners[2] = make_pair(intersections[i]->x, intersections[i]->y);
            idx = i;
        }
    }
    intersections.erase(intersections.begin()+idx);

    dist = maxDistance;
    for (int i = 0; i < intersections.size(); ++i) {
        if (distance(intersections[i]->x - width_, intersections[i]->y - height_) < dist) {
            dist = distance(intersections[i]->x - width_, intersections[i]->y - height_);
            corners[3] = make_pair(intersections[i]->x, intersections[i]->y);
            idx = i;
        }
    }
    intersections.erase(intersections.begin()+idx);

	return corners;
}

CImg<float> CImgUtil::a4_wraping(array<pair<int,int>, 4> corners)
{
    // hardcode size
    int a4_width = 1050, a4_height = 1485;
    CImg<float> a4(a4_width, a4_height, 1, 3);

    // a silly but accurate way to determine portrait or landscape
    bool isPortrait;
    if (distance(corners[0].first - corners[1].first, \
        corners[0].second - corners[1].second) < \
        distance(corners[0].first - corners[2].first, \
        corners[0].second - corners[2].second))
        isPortrait = true;
    else
        isPortrait = false;

    // projection matrix
    float a11, a12, a13, a21, a22, a23, a31, a32, a33;
    int x0, x1, x2, x3, y0, y1, y2, y3;
    if (isPortrait) {
        x0 = corners[0].first;
        x1 = corners[1].first;
        x2 = corners[2].first;
        x3 = corners[3].first;
        y0 = corners[0].second;
        y1 = corners[1].second;
        y2 = corners[2].second;
        y3 = corners[3].second;
    }
    else {
        x0 = corners[2].first;
        x1 = corners[0].first;
        x2 = corners[3].first;
        x3 = corners[1].first;
        y0 = corners[2].second;
        y1 = corners[0].second;
        y2 = corners[3].second;
        y3 = corners[1].second;
    }

    a13 = x0;
    a23 = y0;
    a33 = 1.0f;

    CImg<float> A(2, 2, 1, 1);
    CImg<float> B(1, 2, 1, 1);
    CImg<float> X(1, 2, 1, 1);

    // hardcode the solution
    A.atXY(0, 0) = a4_width*(x3-x1);
    A.atXY(1, 0) = a4_height*(x3-x2);
    A.atXY(0, 1) = a4_width*(y3-y1);
    A.atXY(1, 1) = a4_height*(y3-y2);
    B.atXY(0, 0) = x1-x3+x2-x0;
    B.atXY(0, 1) = y1-y3+y2-y0;
    X = B.get_solve(A);
    a31 = X.atXY(0, 0);
    a32 = X.atXY(0, 1);

    a11 = (x1*(a4_width*a31+1)-a13)/a4_width;
    a21 = (y1*(a4_width*a31+1)-a23)/a4_width;
    a12 = (x2*(a4_height*a32+1)-a13)/a4_height;
    a22 = (y2*(a4_height*a32+1)-a23)/a4_height;

    // inverse projection
    cimg_forXY(a4, u, v) {
        static float raw_x, raw_y;
        raw_x = (a11*u + a12*v + a13) / (a31*u + a32*v + a33);
        raw_x = (raw_x-floor(raw_x)) > 0.5f ? ceil(raw_x):floor(raw_x);
        raw_y = (a21*u + a22*v + a23) / (a31*u + a32*v + a33);
        raw_y = (raw_y-floor(raw_y)) > 0.5f ? ceil(raw_y):floor(raw_y);
        // to keep hard edges
        for (int channel = 0; channel < 3; ++channel) {
            a4.atXY(u, v, channel) = img_rgb_.atXY(int(raw_x), int(raw_y), channel);
        }
    }

    return a4;
}

array<pair<int,int>, 4> CImgUtil::four_corners(bool display, bool save)
{
    return this->four_corners(5, 15, 650, 200, 1, true, 25, display, save);
}

CImg<float> CImgUtil::get_img_rgb()
{
	return img_rgb_;
}

int CImgUtil::get_width()
{
	return width_;
}

int CImgUtil::get_height()
{
	return height_;
}

void CImgUtil::histeq_one(const CImg<float>& src, CImg<float>& img, int dim)
{
	static float bin[256];

	memset(bin, 0, sizeof(float) * 256);

	float factor = 1.0 / (width_ * height_);
	cimg_forXY(src, x, y)
	{
		bin[(int)src(x, y, dim)] += factor;
	}

	for (int i = 1; i < 256; i++)
		bin[i] += bin[i - 1];

	for (int i = 0; i < 256; i++)
		bin[i] = 255.0f * bin[i] + 0.5f;

	cimg_forXY(src, x, y)
	{
		img(x, y, dim) = bin[(int)src(x, y, dim)];
	}
}

void CImgUtil::histeq_intensity(const CImg<float>& src, CImg<float>& img)
{
	static float bin[256];
	int dim = 2; // intensity dimension

	memset(bin, 0, sizeof(float) * 256);

	float factor = 1.0 / (width_ * height_);
	cimg_forXY(src, x, y)
	{
		bin[(int)ceil(src(x, y, dim)*255)] += factor;
	}

	for (int i = 1; i < 256; i++)
		bin[i] += bin[i - 1];

	for (int i = 0; i < 256; i++)
		bin[i] = 255.0f * bin[i] + 0.5f;

	cimg_forXY(src, x, y)
	{
		img(x, y, dim) = bin[(int)ceil(src(x, y, dim)*255)]/255.0;
	}
}

// Calculate the distance
inline float CImgUtil::distance(float x, float y) {
    return sqrt(x*x + y*y);
}

// Polar coordinate intersection at x
inline const int CImgUtil::CrossX(int theta, int distance, int x) {
    float angle = (float)theta*PI / 180.0;
    float m = -cos(angle) / sin(angle);
    float b = (float)distance / sin(angle);
    return m*x + b;
}

// Polar coordinate intersection at y
inline const int CImgUtil::CrossY(int theta, int distance, int y) {
    float angle = (float)theta*PI / 180.0;
    float m = -cos(angle) / sin(angle);
    float b = (float)distance / sin(angle);
    return ((float)(y - b) / m);
}

CImg<float> CImgUtil::non_maximum_suppression(const CImg<float>& response,
											  const array<float, 3>& threshold) {
	CImg<float> result(response);
    result.fill(0);
    for (int i = 1; i < result.width() - 1; i++) {
        for (int j = 1; j < result.height() - 1; j++) {
            for (int k = 0; k < 3; k++) {
                if (response(i, j, k) > threshold[k] && \
                	response(i, j, k) > response(i - 1, j - 1, k) && \
                	response(i, j, k) > response(i - 1, j, k) && \
                	response(i, j, k) > response(i - 1, j + 1, k) && \
                	response(i, j, k) > response(i, j - 1, k) && \
                    response(i, j, k) > response(i, j + 1, k) && \
                    response(i, j, k) > response(i + 1, j - 1, k) && \
                    response(i, j, k) > response(i + 1, j, k) && \
                    response(i, j, k) > response(i + 1, j + 1, k)) {
                    result(i, j, k) = 1;
                }
            }
        }
    }
    return result;
}

CImg<float> CImgUtil::gaussian_lowpass(const CImg<float>& Ixx, \
									   const CImg<float>& Ixy, \
									   const CImg<float>& Iyy, \
									   const float alpha, \
									   array<float, 3>& threshold) {
	CImg<float> response(width_, height_, 1, 3);
    float maxR[3] = { INT_MIN };
    int gauss[3][3] = { {1,2,1},{2,4,2},{1,2,1} }; // assume sigma = 1.0
    for (int i = 1; i < width_-1; i++) {
        for (int j = 1; j < height_-1; j++) {
            for (int k = 0; k < 3; k++) {
                float A = 0.0, B = 0.0, C = 0.0;
                for (int row = 0; row < 3; row++) {
                    for (int col = 0; col < 3; col++) {
                        A += Ixx(row + i - 1, col + j - 1, k)*gauss[row][col];
                        B += Ixy(row + i - 1, col + j - 1, k)*gauss[row][col];
                        C += Iyy(row + i - 1, col + j - 1, k)*gauss[row][col];
                    }
                }
                A /= 16.0; B /= 16.0; C /= 16.0;

                // harris response
                response(i, j, k) = (A*C - B*B) - alpha*(A + C)*(A + C);
                if (maxR[k] < response(i, j, k)) {
                    maxR[k] = response(i, j, k);
                    threshold[k] = 0.01*maxR[k];
                }
            }
        }
    }
    return response;
}