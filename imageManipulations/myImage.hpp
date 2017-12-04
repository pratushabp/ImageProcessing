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


using namespace std;

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
	

	int ditherI2[4] = {1,2,3,0};
	int ditherI4[16] = { 5, 9,6,10,13,1,14,2,7,11,4,8,15,3,12, 0};
	int ditherI8[64] = { 21,37,25,41,22,38,26,42, 53,5,57,9,54,6,58,10 , 29,45,17,33,30,46,18,34 ,61,13,49,1,62,14,50,2 ,23,39,27,43,20,36,24,40 ,55,7,59,11,52,4,56,8 , 31,47,19,35,28,44,16,32 ,63,15,51,3,60,12,48,0  };
	int ditherA4[16] = { 14,10,11,15,9,3,0,4,8,2,1,5,13,7,6,12 };
	//int dither4[16] = { 0, 8, 2, 10, 12, 4, 11, 6, 3, 11, 1, 9, 15, 7 ,13 ,5};

	std::vector<unsigned char>imageData = std::vector<unsigned char>(imageHeight * imageBreadth * bytesPerPixel);
	std::vector<unsigned char>imageDataBinary = std::vector<unsigned char>(imageHeight * imageBreadth);
	std::vector<unsigned char>imageDataOut = std::vector<unsigned char>(imageHeight * imageBreadth);
	std::vector <float>  normalizedData = std::vector<float>(imageHeight * imageBreadth * bytesPerPixel);
	std::vector<unsigned char>fitData = std::vector<unsigned char>(imageHeight * imageBreadth * bytesPerPixel);
	std::vector<unsigned char>outImageData = std::vector<unsigned char>(imageHeight * imageBreadth * bytesPerPixel);
	std::vector<unsigned char>MI = std::vector<unsigned char>(imageHeight * imageBreadth * bytesPerPixel);
	std::vector<unsigned char>MII = std::vector<unsigned char>(imageHeight * imageBreadth * bytesPerPixel);
	std::vector<unsigned char>pieceImage = std::vector<unsigned char>(pieceHeight * pieceBreadth * bytesPerPixel);
	
	//Shrinking 
	unsigned char check[9];
	int bondNum = 0;
	int match;
	int count = 0;
	//[FIXME] Make default constructor
	myImage(int height, int width, int bytesPerPixel, char* inFile, char* outFile) :imageHeight(height), imageBreadth(width), bytesPerPixel(bytesPerPixel), inFile(inFile), outFile(outFile){};
	myImage(int height, int width, int bytesPerPixel, char* inFile, char* outFile, int pieceheight, int piecewidth,  char* inPieceFile) :imageHeight(height), imageBreadth(width), bytesPerPixel(bytesPerPixel), inFile(inFile), outFile(outFile), pieceHeight(pieceheight), pieceBreadth(piecewidth), inPieceFile(inPieceFile){};
	~myImage() {};
	void myImage::readImage();
	void myImage::readPieceImage(int height, int width, int bytes, char* inPieceFile);
	void myImage::writeImage();
	void myImage::writeImage(unsigned char toWrite[], char* inFile, int height, int breadth, int bytes);
	void myImage::nomalize();
	void myImage::convertToBinary(unsigned char imageData[], unsigned char imageDataBinary[]);
	void myImage::binarize(unsigned char imageDataBinary[], unsigned char imageDataOut[]);
	void myImage::medianFilter(unsigned char imageDataBinary[], unsigned char imageDataOut[], int N);
	void myImage::doShrinking(unsigned char imageDataOut[]);
	void myImage::doThinning(unsigned char imageDataOut[]);
	void myImage::doSkel(unsigned char imageDataOut[]);
	void myImage::convertFromDisplay(unsigned char imageDataOut[]);
	int myImage::convertToDisplay(unsigned char imageDataOut[]);
	void myImage::puzzleHillary();
	void myImage::puzzleTrump();
	void myImage::project();
	void myImage::ditheringMatix();
	void myImage::quantizeFourLevels();
	void myImage::errorDiffuseFS();
	void myImage::errorDiffusionJJN();
	void myImage::errorDiffusionStucki();
	void myImage::writePieceImage();
	void myImage::writeImageBinary();
	void myImage::fillHole();
	void myImage::shrinking();
	void myImage::thinning();
	void myImage::skeletonizing();
	void myImage::process();
	unsigned char myImage::hitOrMissStageIShrink(unsigned char arr[], int len);
	unsigned char myImage::hitOrMissStageIThin(unsigned char arr[], int len);
	unsigned char hitOrMissStageISkel(unsigned char arr[], int len);
	unsigned char myImage::hitOrMissStageII(unsigned char arr[], int len);
};
