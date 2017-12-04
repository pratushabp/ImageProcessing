#ifndef _MYIMAGE_H_
#define _MYIMAGE_H_
#include <stdio.h>
#include <iostream>
#include <conio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <new>
#include <iterator>
#include <algorithm>
#include <math.h>
#include <iterator>
#include <algorithm>
#include <math.h>
#include <map>
#include <random>
#include <string.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2\xfeatures2d\nonfree.hpp"
#include "opencv2\xfeatures2d.hpp"
//#include "kmeans.c"

using namespace std;
using namespace cv;

class myImage
{
public:
	int imageHeight;
	int imageBreadth;
	int pieceHeight;
	int pieceBreadth;
	int bytesPerPixel;
	int row_it;
	int col_it;
	int pixel_it;
	int valFrom;
	int valFromR;
	int valFromG;
	int valFromB;
	int valTo;
	char* inFile = 0;
	char* outFile = 0;
	char* inPieceFile = 0;
	int imageExtend = 2;
	int N = 2;

	Mat imageDataMat;

	std::vector<unsigned char>imageData = std::vector<unsigned char>(imageHeight * imageBreadth * bytesPerPixel);
	std::vector<unsigned char>subtractedImage = std::vector<unsigned char>(imageHeight * imageBreadth * bytesPerPixel);
	std::vector<unsigned char>boundaryImage = std::vector<unsigned char>((imageHeight + (2* imageExtend)) * (imageBreadth + (2 * imageExtend)) * bytesPerPixel);
	std::vector<unsigned char>imageDataBinary = std::vector<unsigned char>(imageHeight * imageBreadth);
	std::vector<unsigned char>imageDataOut = std::vector<unsigned char>(imageHeight * imageBreadth);
	std::vector <float>  normalizedData = std::vector<float>(imageHeight * imageBreadth * bytesPerPixel);
	

	//[FIXME] Make default constructor
	myImage(int height, int width, int bytesPerPixel, char* inFile) :imageHeight(height), imageBreadth(width), bytesPerPixel(bytesPerPixel), inFile(inFile) {};
	//myImage(int height, int width, int bytesPerPixel, char* inFile, char* outFile, int pieceheight, int piecewidth, char* inPieceFile) :imageHeight(height), imageBreadth(width), bytesPerPixel(bytesPerPixel), inFile(inFile), outFile(outFile), pieceHeight(pieceheight), pieceBreadth(piecewidth), inPieceFile(inPieceFile) {};
	~myImage() {};
	void myImage::readImage(char inFile[], int imageHeight, int imageBreadth, int bytesPerPixel);
	void myImage::writeImage(unsigned char toWrite[], char* outFile, int height, int breadth, int bytes);
	vector<double> myImage::textClassify();
	vector<vector<unsigned char>> myImage::textSegment();
	vector<vector<unsigned char>> myImage::adtextSegment();
	void myImage::subtractMean();
	void myImage::extendImage();
	void cannyEdge();
	vector<unsigned char>  myImage::applyFilter(double filter[]);
};
#endif