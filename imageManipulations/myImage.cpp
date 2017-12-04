
#include "myImage.hpp"
#include <math.h>

#define M_PI 3.141592653589793

ofstream rankOut("rank.txt");

struct coord {
	int x = 0;
	int y = 0;
};

void myImage::readImage()
{
	FILE* imageRead = 0;
	if (fopen_s(&imageRead,inFile, "rb"))
	{
		cout << " cannot read " << inFile << endl;
	}
	fread(&imageData[0], sizeof(unsigned char), imageHeight * imageBreadth * bytesPerPixel, imageRead);
	



}

void myImage::readPieceImage(int height, int breadth, int bytes, char * inFile)
{

	int count = 0;
	FILE* pieceRead = 0;
	ofstream out("check.txt");
	if (fopen_s(&pieceRead, inPieceFile, "rb"))
	{
		cout << " cannot read " << inFile << endl;
	}
	fread(&pieceImage[0], sizeof(unsigned char), pieceHeight * pieceBreadth * bytesPerPixel, pieceRead);
	
	

}



void myImage::writeImage()
{
	FILE* outImage = 0;
	if (fopen_s(&outImage, outFile, "wb"))
	{
		cout << "Cannot open input file" << endl;
		exit(0);
	}
	fwrite(&outImageData[0], sizeof(unsigned char), imageHeight * imageBreadth * bytesPerPixel, outImage);
}

void myImage::writeImage(unsigned char toWrite[], char* inFile, int height, int breadth, int bytes)
{

	FILE* outImage = 0;
	if (fopen_s(&outImage, inFile, "wb"))
	{
		cout << "Cannot open input file" << endl;
		exit(0);
	}
	fwrite(toWrite, sizeof(unsigned char), height * breadth * bytes, outImage);
}


void myImage::writePieceImage()
{
	FILE* outImage = 0;
	if (fopen_s(&outImage, "piece_out.raw", "wb"))
	{
		cout << "Cannot open input file" << endl;
		exit(0);
	}
	fwrite(&pieceImage[0], sizeof(unsigned char), pieceHeight * pieceBreadth * bytesPerPixel, outImage);
}


void myImage::writeImageBinary()
{
	
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			valFrom = row_it * imageBreadth + col_it;
			imageDataOut[valFrom] *= 255;
			
		}
	}
	
	FILE* outImage = 0;
	if (fopen_s(&outImage, outFile, "wb"))
	{
		cout << "Cannot open input file" << endl;
		exit(0);
	}
	fwrite(&imageDataOut[0], sizeof(unsigned char), imageHeight * imageBreadth , outImage);
}






void myImage::nomalize()
{
	
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
			{
				valFrom = (row_it * imageBreadth + col_it) * bytesPerPixel + pixel_it; // grey scale and color normalisation
				valTo = valFrom;
				normalizedData[valTo] = double(imageData[valFrom] )/ 255;
			}
		}
	}
}

void myImage::puzzleHillary()
{
	int pieceHeight = 500; int pieceWidth = 500;

	coord topCorner;
	coord botCorner;
	coord leftCorner;
	coord rightCorner;

	int rowMin = 512;
	int rowMax = 0;
	int colMin = 512;
	int colMax = 0;
	int cutHeight;
	int cutBreadth;

	//Covert target image to grey
	std::vector<unsigned char>pieceImageBinary = std::vector<unsigned char>(pieceHeight * pieceBreadth);
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			imageDataBinary[row_it * imageBreadth + col_it] = (imageData[(row_it * imageBreadth + col_it) * 3] + imageData[(row_it * imageBreadth + col_it) * 3 + 1] + imageData[(row_it * imageBreadth + col_it) * 3 + 2]) / 3;
		}
	}

	//Covert piece to grey
	for (row_it = 0; row_it < pieceHeight; row_it++)
	{
		for (col_it = 0; col_it < pieceBreadth; col_it++)
		{

			pieceImageBinary[row_it * pieceBreadth + col_it] = (pieceImage[(row_it * pieceBreadth + col_it) * 3] + pieceImage[(row_it * pieceBreadth + col_it) * 3 + 1] + pieceImage[(row_it * pieceBreadth + col_it) * 3 + 2]) / 3;
		}
	}



	for (int i = 0; i < imageHeight; i++)
	{
		for (int j = 0; j < imageBreadth; j++)
		{
			//out << j << " " << i << endl;
			unsigned char pixel_value = imageDataBinary[(i * imageBreadth + j)];
			if (pixel_value == 255)
			{
				if (i < rowMin)
				{
					rowMin = i;
					topCorner.y = i;
					topCorner.x = j;
				}
				if (i > rowMax)
				{
					rowMax = i;
					botCorner.y = i;
					botCorner.x = j;
				}
				if (j < colMin) {
					colMin = j;
					leftCorner.y = i;
					leftCorner.x = j;
				}
				if (j > colMax) {
					colMax = j;
					rightCorner.y = i;
					rightCorner.x = j;
				}
			}
		}
	}

	cutHeight = botCorner.y - topCorner.y;
	cutBreadth = rightCorner.x - leftCorner.x;



	std::cout << "topLeft " << colMin << " " << rowMin << std::endl;
	std::cout << "botLeft " << colMin << " " << rowMax << std::endl;
	std::cout << "topRight " << colMax << " " << rowMin << std::endl;
	std::cout << " botRight" << colMax << " " << rowMax << std::endl << std::endl;


	//writeImage(&imageDataBinary[0], 512, 512);


	//finding the corners of the piece image
	rowMin = pieceHeight;
	coord pieceTopCorner;
	rowMax = 0;
	coord pieceBotCorner;
	colMin = pieceWidth;
	coord pieceLeftCorner;
	colMax = 0;
	coord pieceRightCorner;
	int temp;
	ofstream out("isiteratingright.txt");
	for (int i = 0; i < pieceHeight / 2; i++)
	{
		for (int j = 0; j < pieceWidth; j++)
		{
			//out << j << " " << i << endl
			int pixel_value = pieceImageBinary[(i * pieceBreadth + j)];
			if (pixel_value != 255)
			{
				if (i < rowMin)
				{
					rowMin = i;
					pieceTopCorner.y = i;
					pieceTopCorner.x = j;
				}
				if (i > rowMax)
				{
					rowMax = i;
					pieceBotCorner.y = i;
					pieceBotCorner.x = j;
				}
				if (j < colMin) {
					colMin = j;
					pieceLeftCorner.y = i;
					pieceLeftCorner.x = j;
				}
				if (j > colMax) {
					colMax = j;
					pieceRightCorner.y = i;
					pieceRightCorner.x = j;
				}
			}
		}
	}


	std::cout << "leftCorner of piece:" << pieceLeftCorner.x << " " << pieceLeftCorner.y << std::endl;
	std::cout << "topCorner of piece:" << pieceTopCorner.x << " " << pieceTopCorner.y << std::endl;
	std::cout << "botCorner of piece:" << pieceBotCorner.x << " " << pieceBotCorner.y << std::endl;
	std::cout << "rightCorner of piece:" << pieceRightCorner.x << " " << pieceRightCorner.y << std::endl << std::endl;

	/*std::cout << "x_min:" << x_min << std::endl;
	std::cout << "x_max:" << x_max << std::endl;
	std::cout << "y_min:" << y_min << std::endl;
	std::cout << "y_max:" << y_max << std::endl;*/

	cutHeight = pieceBotCorner.y - pieceTopCorner.y;
	cutBreadth = pieceRightCorner.x - pieceLeftCorner.x;

	std::vector<unsigned char>debugImage = std::vector<unsigned char>(cutHeight * cutBreadth * bytesPerPixel);
	int debugCount = 0;
	for (row_it = pieceTopCorner.y; row_it < pieceBotCorner.y; row_it++)
	{
		for (col_it = pieceLeftCorner.x; col_it < pieceRightCorner.x; col_it++)
		{
			for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
			{

				debugImage[debugCount++] = pieceImage[(row_it * pieceBreadth + col_it) * bytesPerPixel + pixel_it];
			}
		}
	}


	/*double H[3][3] = {
	0.6673,
	0.1491,
	100.9746,
	- 0.1524,
	0.6083,
	115.3751,
	0.0001,
	- 0.0002,
	1};*/

	double H[3][3] = { 0.667319093994468	,0.149053569311279	,100.974620419000,
		-0.152425172540027	,0.608313239330649,	115.375126260629,
		0.000124941009088452 ,-0.000156344527505222,	1 };

	double HInv[3][3] = { 0.667319093994468 ,-0.152425172540027,	0.000124941009088452,
		0.149053569311279,	0.608313239330649 ,-0.000156344527505222,
		100.974620419000,	115.375126260629	,1 };


	//double H[3][3] = { 0.667319093994468 , -0.152425172540027	,0.000124941009088448,		0.149053569311279,	0.608313239330649, -0.000156344527505218,		100.974620419000	,115.375126260629	,1 };

	//double H[3][3] = { 1.46118979901988	,0.389214807447154 ,- 0.000121710922790932,
	//- 0.384549419153414	,1.52733035241520,	0.000286835734766136,
	//- 103.175647527870 ,- 215.516749695026,	1 };


	double xIndNew = 0, yIndNew = 0, xIndNext, yIndNext, xFloor, yFloor;
	int rowNext, colNext;
	double xCart, yCart, xCartNew, yCartNew = 0;
	unsigned char value, pixel_value;
	debugCount = 0;
	for (row_it = rowMin; row_it < rowMax; row_it++)
	{
		for (col_it = colMin; col_it < colMax; col_it++)
		{
			pixel_value = pieceImageBinary[(row_it * pieceBreadth + col_it)];
			//pixel_value = (pieceImage[(row_it * imageBreadth + col_it) * bytesPerPixel + 0] + pieceImage[(row_it * imageBreadth + col_it) * bytesPerPixel + 1] + pieceImage[(row_it * imageBreadth + col_it) * bytesPerPixel + 2]) / 3;
			if (pixel_value != 255)
			{
				xIndNew = (H[0][0] * col_it + H[0][1] * row_it + H[0][2]) / (H[2][0] * col_it + H[2][1] * row_it + 1);
				yIndNew = (H[1][0] * col_it + H[1][1] * row_it + H[1][2]) / (H[2][0] * col_it + H[2][1] * row_it + 1);
				xFloor = floor(xIndNew);
				yFloor = floor(yIndNew);
				xIndNext = xFloor + 1;
				if (xIndNext > cutBreadth)
				{
					xIndNext = cutBreadth - 1;
				}
				yIndNext = yFloor + 1;
				if (yIndNext > cutHeight)
				{
					yIndNext = cutHeight - 1;
				}
				//if (xFloor > x_min && xFloor < x_max && yFloor > y_min && yFloor < y_max)
				{
					for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
					{
						value = pieceImage[(row_it * pieceBreadth + col_it) * bytesPerPixel + pixel_it];
						out << "Indx " << round(xIndNew) << " " << round(yIndNew) << "  " << col_it << "  " << row_it << "  " << value << endl;
						imageData[(yFloor * imageBreadth + xFloor) * bytesPerPixel + pixel_it] = value;
					}
				}
			}
		}
	}

	row_it = topCorner.y;
	{
		for (int i = row_it - 2; i < row_it + 2; i++)
		{
			for (int col_it = topCorner.x; col_it < rightCorner.x; col_it++)
			{
				for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
				{
					imageData[(i * imageBreadth + col_it) * bytesPerPixel + pixel_it] = ((imageData[((i + 5) * imageBreadth + col_it) * bytesPerPixel + pixel_it] + imageData[((i - 5) * imageBreadth + col_it)* bytesPerPixel + pixel_it]) / 2);
				}
			}
		}
	}

	row_it = botCorner.y;
	{
		for (int i = row_it - 2; i < row_it + 2; i++)
		{
			for (int col_it = topCorner.x; col_it < rightCorner.x; col_it++)
			{
				for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
				{
					imageData[(i * imageBreadth + col_it) * bytesPerPixel + pixel_it] = ((imageData[((i + 4) * imageBreadth + col_it) * bytesPerPixel + pixel_it] + imageData[((i - 4) * imageBreadth + col_it)* bytesPerPixel + pixel_it]) / 2);
				}
			}
		}
	}

	//int col_it = 173;

	for (row_it = 135; row_it < 234; row_it++)
	{
		for (col_it = 173 - 2; col_it < 173 + 2; col_it++) {
			for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
			{
				imageData[(row_it * imageBreadth + col_it) * bytesPerPixel + pixel_it] = (imageData[(row_it * imageBreadth + (col_it - 4))* bytesPerPixel + pixel_it] + imageData[(row_it * imageBreadth + (col_it + 4))* bytesPerPixel + pixel_it]) / 2;

			}
		}
	}


	for (row_it = 135; row_it <= 234; row_it++)
	{
		for (col_it = 272 - 5; col_it < 272 + 5; col_it++) {
			for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
			{
				imageData[(row_it * imageBreadth + col_it) * bytesPerPixel + pixel_it] = (imageData[(row_it * imageBreadth + (col_it - 9))* bytesPerPixel + pixel_it] + imageData[(row_it * imageBreadth + (col_it + 9))* bytesPerPixel + pixel_it]) / 2;

			}
		}
	}


	writeImage(&imageData[0],"cut_hillary.raw",512, 512, 3);


	/*
	std::vector<unsigned char>debugImage = std::vector<unsigned char>(cutHeight * cutBreadth * bytesPerPixel);
	int debugCount = 0;
	for (row_it = pieceTopLeft.y; row_it < pieceBotCorner.y; row_it++)
	{
	for (col_it = pieceLeftCorner.x; col_it < pieceRightCorner.x; col_it++)
	{
	for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
	{
	debugImage[debugCount++] = pieceImage[(row_it * pieceBreadth + col_it) * bytesPerPixel + pixel_it];
	}
	}
	}
	int count = 0;
	double value;
	//writeImage(&debugImage[0], cutHeight, cutBreadth);
	//out << y_max - y_min << " " << x_max - x_min << endl;
	std::vector<unsigned char>fitData = std::vector<unsigned char>(pieceHeight * pieceBreadth * bytesPerPixel);
	//rotation of image
	//angle value found by trial and error
	double angle = 14.5 / 180.0 * M_PI;
	//double angle = atan((topCorner.x - leftCorner.x) / (topCorner.y - leftCorner.y));
	for (int y = y_min; y < y_max; y++) {
	for (int x = x_min; x < x_max; x++) {
	for (int k = 0; k < bytesPerPixel; k++) {
	//std::cout << x << " " << y << std::endl;
	double x_new = x - x_max / 2;
	double y_new = y - y_max / 2;
	x_new = cos(angle) * x_new - sin(angle) * y_new;
	y_new = sin(angle) * x_new + cos(angle) * y_new;
	x_new = x_new + x_max / 2;
	y_new = y_new + y_max / 2;
	//std::cout << x_new << " " << y_new << std::endl;
	x_new = (x_new);
	y_new = (y_new);
	out << count++ << endl;
	if (x_new > 0 && y_new > 0)
	{
	fitData[(y * pieceBreadth + x) * bytesPerPixel + k] = pieceImage[(int(y_new) * pieceBreadth + int(x_new)) * bytesPerPixel + k];

	}
	}
	}
	}

	//writeImage(&fitData[0], 500, 500);
	//write(imageMod, "trumpOut.raw", 500, 500);


	/*
	y_min = 500;
	y_max = 0;
	x_min = 500;
	x_max = 0;

	for (int i = 0; i < pHeight / 2; i++)
	{
	for (int j = 0; j < pWidth; j++)
	{
	//out << j << " " << i << endl;
	for (int k = 0; k < bytesPerPixel; k++)
	{
	temp = (i * pieceBreadth + j) * bytesPerPixel + k;

	int pixel_value = fitData[(i * pieceBreadth + j) * bytesPerPixel + k];
	if (pixel_value != 255 && pixel_value != 0)
	{
	if (i < y_min)
	{
	y_min = i;
	topCorner.y = i;
	topCorner.x = j;
	}
	if (i > y_max)
	{
	y_max = i;
	botCorner.y = i;
	botCorner.x = j;
	}
	if (j < x_min) {
	x_min = j;
	leftCorner.y = i;
	leftCorner.x = j;
	}
	if (j > x_max) {
	x_max = j;
	rightCorner.y = i;
	rightCorner.x = j;
	}
	}
	}
	}
	}
	cutHeight = botCorner.y - topCorner.y;
	cutBreadth = rightCorner.x - leftCorner.x;
	std::vector<unsigned char>debugRotated = std::vector<unsigned char>(cutHeight * cutBreadth * bytesPerPixel);

	debugCount = 0;
	for (row_it = topCorner.y; row_it < botCorner.y; row_it++)
	{
	for (col_it = leftCorner.x; col_it < rightCorner.x; col_it++)
	{
	for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
	{
	debugRotated[debugCount++] = fitData[(row_it * pieceBreadth + col_it) * bytesPerPixel + pixel_it];
	}
	}
	}



	*/

	//	writeImage(&debugRotated[0], cutHeight, cutBreadth);

	/*std::vector<unsigned char>downScaleData = std::vector<unsigned char>(imageHeight * imageBreadth * bytesPerPixel);
	//Down-scaling of image piece

	for (int y = y_min; y < y_max; y++) {
	for (int x = x_min; x < x_max; x++) {
	for (int k = 0; k < bytesPerPixel; k++) {
	if (fitData[(y * imageBreadth + x) * bytesPerPixel + k] != 255) {
	//std::cout << x << " " << y << std::endl;
	double x_new = x * 0.69;
	double y_new = y * 0.6666;
	//std::cout << x_new << " " << y_new << std::endl;
	downScaleData[(round(y_new) * imageBreadth + round(x_new)) * bytesPerPixel + k] = fitData[(y * pieceBreadth + x) * bytesPerPixel + k];
	}
	}
	}
	}

	//writeImage(&downScaleData[0], 512, 512);
	//	ofstream out("whereiszero.txt");
	//finding the corners of the missing piece in the original image




	y_min = imageHeight;
	y_max = 0;
	x_min = imageBreadth;
	x_max = 0;

	for (int i = 0; i < imageHeight; i++) {
	for (int j = 0; j < imageBreadth; j++) {
	out << i * imageBreadth + j << endl;
	int pixel_value = (imageData[i * imageBreadth + j + 0] + imageData[i * imageBreadth + j + 1] + imageData[i * imageBreadth + j + 2]) / 3;
	//	out << pixel_value << endl;
	if (pixel_value == 255)
	{
	count = 0;
	for (int p = 0; p < 100; p++)
	{
	if (imageData[i * imageBreadth + j + p] == 255)
	++count;
	}


	if (count > 90)
	{
	topCorner.x = j;
	topCorner.y = i;
	break;
	}
	/*if (i < y_min) {
	y_min = i;
	topCorner.y = i;
	topCorner.x = j;
	}
	else if (i >= y_max) {
	y_max = i;
	botCorner.y = i;
	botCorner.x = j;
	}
	if (j <= x_min) {
	x_min = j;
	leftCorner.y = i;
	leftCorner.x = j;
	}
	else if (j > x_max) {
	x_max = j;
	rightCorner.y = i;
	rightCorner.x = j;*//*
	}

	}
	}



	//	writeImage(&imageData[0], 512, 512);
	/*
	std::cout << "leftCorner:" << leftCorner.x << " " << leftCorner.y << std::endl;
	std::cout << "topCorner:" << topCorner.x << " " << topCorner.y << std::endl;
	std::cout << "botCorner:" << botCorner.x << " " << botCorner.y << std::endl;
	std::cout << "rightCorner:" << rightCorner.x << " " << rightCorner.y << std::endl << std::endl;


	//Replacing the missing piece with the transformed piece
	for (int i = y_min; i < y_max + 1; i++) {
	for (int j = x_min; j < x_max + 1; j++) {
	for (int k = 0; k < bytesPerPixel; k++) {
	int i_ind = 43 + i - y_min;
	int j_ind = 57 + j - x_min;
	imageData[(i * imageBreadth + j) * bytesPerPixel + k] = fitData[(i_ind * pieceBreadth + j_ind ) * bytesPerPixel +  k];
	}
	}
	}

	for (int i = y_min; i < y_max + 1; i++) {
	for (int j = x_min; j < x_max + 1; j++) {
	for (int k = 0; k < bytesPerPixel; k++) {
	int i_ind = 43 + i - y_min;
	int j_ind = 57 + j - x_min;
	outImageData[(i * imageBreadth + j) * bytesPerPixel + k] = imageData[(i * imageBreadth + j) * bytesPerPixel + k];
	}
	}
	}
	writeImage(&outImageData[0], 512, 512);*/
}

void myImage::puzzleTrump()

{
	int pieceHeight = 500; int pieceWidth = 500;

	coord topCorner;
	coord botCorner;
	coord leftCorner;
	coord rightCorner;

	int rowMin = 512;
	int rowMax = 0;
	int colMin = 512;
	int colMax = 0;
	int cutHeight;
	int cutBreadth;

	//Covert target image to grey
	std::vector<unsigned char>pieceImageBinary = std::vector<unsigned char>(pieceHeight * pieceBreadth);
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			imageDataBinary[row_it * imageBreadth + col_it] = (imageData[(row_it * imageBreadth + col_it) * 3] + imageData[(row_it * imageBreadth + col_it) * 3 + 1] + imageData[(row_it * imageBreadth + col_it) * 3 + 2]) / 3;
		}
	}

	//Covert piece to grey
	for (row_it = 0; row_it < pieceHeight; row_it++)
	{
		for (col_it = 0; col_it < pieceBreadth; col_it++)
		{

			pieceImageBinary[row_it * pieceBreadth + col_it] = (pieceImage[(row_it * pieceBreadth + col_it) * 3] + pieceImage[(row_it * pieceBreadth + col_it) * 3 + 1] + pieceImage[(row_it * pieceBreadth + col_it) * 3 + 2]) / 3;
		}
	}



	for (int i = 0; i < imageHeight; i++)
	{
		for (int j = 0; j < imageBreadth; j++)
		{
			//out << j << " " << i << endl;
			unsigned char pixel_value = imageDataBinary[(i * imageBreadth + j)];
			if (pixel_value == 255)
			{
				if (i < rowMin)
				{
					rowMin = i;
					topCorner.y = i;
					topCorner.x = j;
				}
				if (i > rowMax)
				{
					rowMax = i;
					botCorner.y = i;
					botCorner.x = j;
				}
				if (j < colMin) {
					colMin = j;
					leftCorner.y = i;
					leftCorner.x = j;
				}
				if (j > colMax) {
					colMax = j;
					rightCorner.y = i;
					rightCorner.x = j;
				}
			}
		}
	}

	cutHeight = botCorner.y - topCorner.y;
	cutBreadth = rightCorner.x - leftCorner.x;



	std::cout << "topLeft " << colMin << " " << rowMin << std::endl;
	std::cout << "botLeft " << colMin << " " << rowMax << std::endl;
	std::cout << "topRight " << colMax << " " << rowMin << std::endl;
	std::cout << " botRight" << colMax << " " << rowMax << std::endl << std::endl;


	//writeImage(&imageDataBinary[0], 512, 512);


	//finding the corners of the piece image
	rowMin = pieceHeight;
	coord pieceTopCorner;
	rowMax = 0;
	coord pieceBotCorner;
	colMin = pieceWidth;
	coord pieceLeftCorner;
	colMax = 0;
	coord pieceRightCorner;
	int temp;
	ofstream out("isiteratingright.txt");
	for (int i = pieceHeight / 2; i < pieceHeight; i++)
	{
		for (int j = 0; j < pieceWidth; j++)
		{
			//out << j << " " << i << endl
			int pixel_value = pieceImageBinary[(i * pieceBreadth + j)];
			if (pixel_value != 255)
			{
				if (i < rowMin)
				{
					rowMin = i;
					pieceTopCorner.y = i;
					pieceTopCorner.x = j;
				}
				if (i > rowMax)
				{
					rowMax = i;
					pieceBotCorner.y = i;
					pieceBotCorner.x = j;
				}
				if (j < colMin) {
					colMin = j;
					pieceLeftCorner.y = i;
					pieceLeftCorner.x = j;
				}
				if (j > colMax) {
					colMax = j;
					pieceRightCorner.y = i;
					pieceRightCorner.x = j;
				}
			}
		}
	}


	std::cout << "leftCorner of piece:" << pieceLeftCorner.x << " " << pieceLeftCorner.y << std::endl;
	std::cout << "topCorner of piece:" << pieceTopCorner.x << " " << pieceTopCorner.y << std::endl;
	std::cout << "botCorner of piece:" << pieceBotCorner.x << " " << pieceBotCorner.y << std::endl;
	std::cout << "rightCorner of piece:" << pieceRightCorner.x << " " << pieceRightCorner.y << std::endl << std::endl;

	/*std::cout << "x_min:" << x_min << std::endl;
	std::cout << "x_max:" << x_max << std::endl;
	std::cout << "y_min:" << y_min << std::endl;
	std::cout << "y_max:" << y_max << std::endl;*/

	cutHeight = pieceBotCorner.y - pieceTopCorner.y;
	cutBreadth = pieceRightCorner.x - pieceLeftCorner.x;

	std::vector<unsigned char>debugImage = std::vector<unsigned char>(cutHeight * cutBreadth * bytesPerPixel);
	int debugCount = 0;
	for (row_it = pieceTopCorner.y; row_it < pieceBotCorner.y; row_it++)
	{
		for (col_it = pieceLeftCorner.x; col_it < pieceRightCorner.x; col_it++)
		{
			for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
			{

				debugImage[debugCount++] = pieceImage[(row_it * pieceBreadth + col_it) * bytesPerPixel + pixel_it];
			}
		}
	}


	/*double H[3][3] = {
	0.6673,
	0.1491,
	100.9746,
	- 0.1524,
	0.6083,
	115.3751,
	0.0001,
	- 0.0002,
	1};*/

	double H[3][3] = { -0.0495 ,   1.3682,   71.2024,
		-0.9236 ,   0.2406,  558.7460,
		0.0001 ,   0.0011,    0.9810 };

	double HInv[3][3] = { 0.667319093994468 ,-0.152425172540027,	0.000124941009088452,
		0.149053569311279,	0.608313239330649 ,-0.000156344527505222,
		100.974620419000,	115.375126260629	,1 };


	//double H[3][3] = { 0.667319093994468 , -0.152425172540027	,0.000124941009088448,		0.149053569311279,	0.608313239330649, -0.000156344527505218,		100.974620419000	,115.375126260629	,1 };

	//double H[3][3] = { 1.46118979901988	,0.389214807447154 ,- 0.000121710922790932,
	//- 0.384549419153414	,1.52733035241520,	0.000286835734766136,
	//- 103.175647527870 ,- 215.516749695026,	1 };


	double xIndNew = 0, yIndNew = 0, xIndNext, yIndNext, xFloor, yFloor;
	int rowNext, colNext;
	double xCart, yCart, xCartNew, yCartNew = 0;
	unsigned char  pixel_value;
	int value;
	debugCount = 0;
	for (row_it = 236; row_it < 335; row_it++)
	{
		for (col_it = 163; col_it < 262; col_it++)
		{
			pixel_value = pieceImageBinary[(row_it * pieceBreadth + col_it)];
			//pixel_value = (pieceImage[(row_it * imageBreadth + col_it) * bytesPerPixel + 0] + pieceImage[(row_it * imageBreadth + col_it) * bytesPerPixel + 1] + pieceImage[(row_it * imageBreadth + col_it) * bytesPerPixel + 2]) / 3;
			for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
			{
				xIndNew = (H[0][0] * col_it + H[0][1] * row_it + H[0][2]) / (H[2][0] * col_it + H[2][1] * row_it + 1);
				yIndNew = (H[1][0] * col_it + H[1][1] * row_it + H[1][2]) / (H[2][0] * col_it + H[2][1] * row_it + 1);
				xFloor = (xIndNew);
				yFloor = (yIndNew);
				xIndNext = xFloor + 1;
				if (xIndNext > cutBreadth)
				{
					xIndNext = cutBreadth - 1;
				}
				yIndNext = yFloor + 1;
				if (yIndNext > cutHeight)
				{
					yIndNext = cutHeight - 1;
				}
				//if (xFloor > x_min && xFloor < x_max && yFloor > y_min && yFloor < y_max)
				{

					int y_yIndNew = yIndNew;
					int x_xIndNew = xIndNew;

					yIndNew = yIndNew - y_yIndNew;
					xIndNew = xIndNew - x_xIndNew;

					// forward and backward index 
					int r_index = x_xIndNew + 1;
					int b_index = y_yIndNew + 1;

					if (r_index >= pieceBreadth)
						r_index %= imageBreadth;


					if (b_index >= pieceBreadth)
						b_index %= imageBreadth;

					int a, b, c, d;

					a = (1 - yIndNew) * (1 - xIndNew) * pieceImage[(y_yIndNew * pieceBreadth + x_xIndNew) * bytesPerPixel + pixel_it];
					b = (1 - yIndNew) * xIndNew * pieceImage[(y_yIndNew * pieceBreadth + r_index) * bytesPerPixel + pixel_it];
					c = yIndNew * (1 - xIndNew) * pieceImage[(b_index * pieceBreadth + x_xIndNew) * bytesPerPixel + pixel_it];
					d = yIndNew * xIndNew * pieceImage[(b_index * pieceBreadth + r_index) * bytesPerPixel + pixel_it];
					value = (a + b + c + d);
					if (value != (unsigned char)255)
					{
						imageData[(row_it * imageBreadth + col_it) * bytesPerPixel + pixel_it] = value;
					}



					/*value = pieceImage[(row_it * pieceBreadth + col_it) * bytesPerPixel + pixel_it];
					out << "Indx " << round(xIndNew) << " " << round(yIndNew) << "  " << col_it << "  " << row_it << "  " << value << endl;
					imageData[(yFloor * imageBreadth + xFloor) * bytesPerPixel + pixel_it] = value;*/
				}
			}

		}
	}


	//writeImage(&imageData[0], 512, 512, 3);
	row_it = topCorner.y;
	{
		for (int i = row_it; i < row_it + 8; i++)
		{
			for (int col_it = topCorner.x; col_it < rightCorner.x; col_it++)
			{
				for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
				{
					imageData[(i * imageBreadth + col_it) * bytesPerPixel + pixel_it] = ((imageData[((i + 10) * imageBreadth + col_it) * bytesPerPixel + pixel_it] + imageData[((i - 10) * imageBreadth + col_it)* bytesPerPixel + pixel_it]) / 2);
				}
			}
		}
	}

	row_it = botCorner.y;
	{
		for (int i = row_it - 1; i < row_it + 1; i++)
		{
			for (int col_it = topCorner.x; col_it < rightCorner.x; col_it++)
			{
				for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
				{
					imageData[(i * imageBreadth + col_it) * bytesPerPixel + pixel_it] = ((imageData[((i + 2) * imageBreadth + col_it) * bytesPerPixel + pixel_it] + imageData[((i - 2) * imageBreadth + col_it)* bytesPerPixel + pixel_it]) / 2);
				}
			}
		}
	}

	//int col_it = 173;

	for (row_it = 236; row_it < 336; row_it++)
	{
		for (col_it = 262 - 10; col_it < 262 + 3; col_it++) {
			for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
			{
				imageData[(row_it * imageBreadth + col_it) * bytesPerPixel + pixel_it] = (imageData[(row_it * imageBreadth + (col_it - 15))* bytesPerPixel + pixel_it] + imageData[(row_it * imageBreadth + (col_it + 15))* bytesPerPixel + pixel_it]) / 2;

			}
		}
	}

	/*
	for (row_it = 135; row_it <= 234; row_it++)
	{
	for (col_it = 272 - 5; col_it < 272 + 5; col_it++) {
	for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
	{
	imageData[(row_it * imageBreadth + col_it) * bytesPerPixel + pixel_it] = (imageData[(row_it * imageBreadth + (col_it - 9))* bytesPerPixel + pixel_it] + imageData[(row_it * imageBreadth + (col_it + 9))* bytesPerPixel + pixel_it]) / 2;

	}
	}
	}
	*/

	writeImage(&imageData[0], "cut_trump.raw",512, 512, 3);


	/*
	std::vector<unsigned char>debugImage = std::vector<unsigned char>(cutHeight * cutBreadth * bytesPerPixel);
	int debugCount = 0;
	for (row_it = pieceTopLeft.y; row_it < pieceBotCorner.y; row_it++)
	{
	for (col_it = pieceLeftCorner.x; col_it < pieceRightCorner.x; col_it++)
	{
	for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
	{
	debugImage[debugCount++] = pieceImage[(row_it * pieceBreadth + col_it) * bytesPerPixel + pixel_it];
	}
	}
	}
	int count = 0;
	double value;
	//writeImage(&debugImage[0], cutHeight, cutBreadth);
	//out << y_max - y_min << " " << x_max - x_min << endl;
	std::vector<unsigned char>fitData = std::vector<unsigned char>(pieceHeight * pieceBreadth * bytesPerPixel);
	//rotation of image
	//angle value found by trial and error
	double angle = 14.5 / 180.0 * M_PI;
	//double angle = atan((topCorner.x - leftCorner.x) / (topCorner.y - leftCorner.y));
	for (int y = y_min; y < y_max; y++) {
	for (int x = x_min; x < x_max; x++) {
	for (int k = 0; k < bytesPerPixel; k++) {
	//std::cout << x << " " << y << std::endl;
	double x_new = x - x_max / 2;
	double y_new = y - y_max / 2;
	x_new = cos(angle) * x_new - sin(angle) * y_new;
	y_new = sin(angle) * x_new + cos(angle) * y_new;
	x_new = x_new + x_max / 2;
	y_new = y_new + y_max / 2;
	//std::cout << x_new << " " << y_new << std::endl;
	x_new = (x_new);
	y_new = (y_new);
	out << count++ << endl;
	if (x_new > 0 && y_new > 0)
	{
	fitData[(y * pieceBreadth + x) * bytesPerPixel + k] = pieceImage[(int(y_new) * pieceBreadth + int(x_new)) * bytesPerPixel + k];

	}
	}
	}
	}

	//writeImage(&fitData[0], 500, 500);
	//write(imageMod, "trumpOut.raw", 500, 500);


	/*
	y_min = 500;
	y_max = 0;
	x_min = 500;
	x_max = 0;

	for (int i = 0; i < pHeight / 2; i++)
	{
	for (int j = 0; j < pWidth; j++)
	{
	//out << j << " " << i << endl;
	for (int k = 0; k < bytesPerPixel; k++)
	{
	temp = (i * pieceBreadth + j) * bytesPerPixel + k;

	int pixel_value = fitData[(i * pieceBreadth + j) * bytesPerPixel + k];
	if (pixel_value != 255 && pixel_value != 0)
	{
	if (i < y_min)
	{
	y_min = i;
	topCorner.y = i;
	topCorner.x = j;
	}
	if (i > y_max)
	{
	y_max = i;
	botCorner.y = i;
	botCorner.x = j;
	}
	if (j < x_min) {
	x_min = j;
	leftCorner.y = i;
	leftCorner.x = j;
	}
	if (j > x_max) {
	x_max = j;
	rightCorner.y = i;
	rightCorner.x = j;
	}
	}
	}
	}
	}
	cutHeight = botCorner.y - topCorner.y;
	cutBreadth = rightCorner.x - leftCorner.x;
	std::vector<unsigned char>debugRotated = std::vector<unsigned char>(cutHeight * cutBreadth * bytesPerPixel);

	debugCount = 0;
	for (row_it = topCorner.y; row_it < botCorner.y; row_it++)
	{
	for (col_it = leftCorner.x; col_it < rightCorner.x; col_it++)
	{
	for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
	{
	debugRotated[debugCount++] = fitData[(row_it * pieceBreadth + col_it) * bytesPerPixel + pixel_it];
	}
	}
	}



	*/

	//	writeImage(&debugRotated[0], cutHeight, cutBreadth);

	/*std::vector<unsigned char>downScaleData = std::vector<unsigned char>(imageHeight * imageBreadth * bytesPerPixel);
	//Down-scaling of image piece

	for (int y = y_min; y < y_max; y++) {
	for (int x = x_min; x < x_max; x++) {
	for (int k = 0; k < bytesPerPixel; k++) {
	if (fitData[(y * imageBreadth + x) * bytesPerPixel + k] != 255) {
	//std::cout << x << " " << y << std::endl;
	double x_new = x * 0.69;
	double y_new = y * 0.6666;
	//std::cout << x_new << " " << y_new << std::endl;
	downScaleData[(round(y_new) * imageBreadth + round(x_new)) * bytesPerPixel + k] = fitData[(y * pieceBreadth + x) * bytesPerPixel + k];
	}
	}
	}
	}

	//writeImage(&downScaleData[0], 512, 512);
	//	ofstream out("whereiszero.txt");
	//finding the corners of the missing piece in the original image




	y_min = imageHeight;
	y_max = 0;
	x_min = imageBreadth;
	x_max = 0;

	for (int i = 0; i < imageHeight; i++) {
	for (int j = 0; j < imageBreadth; j++) {
	out << i * imageBreadth + j << endl;
	int pixel_value = (imageData[i * imageBreadth + j + 0] + imageData[i * imageBreadth + j + 1] + imageData[i * imageBreadth + j + 2]) / 3;
	//	out << pixel_value << endl;
	if (pixel_value == 255)
	{
	count = 0;
	for (int p = 0; p < 100; p++)
	{
	if (imageData[i * imageBreadth + j + p] == 255)
	++count;
	}


	if (count > 90)
	{
	topCorner.x = j;
	topCorner.y = i;
	break;
	}
	/*if (i < y_min) {
	y_min = i;
	topCorner.y = i;
	topCorner.x = j;
	}
	else if (i >= y_max) {
	y_max = i;
	botCorner.y = i;
	botCorner.x = j;
	}
	if (j <= x_min) {
	x_min = j;
	leftCorner.y = i;
	leftCorner.x = j;
	}
	else if (j > x_max) {
	x_max = j;
	rightCorner.y = i;
	rightCorner.x = j;*//*
	}

	}
	}



	//	writeImage(&imageData[0], 512, 512);
	/*
	std::cout << "leftCorner:" << leftCorner.x << " " << leftCorner.y << std::endl;
	std::cout << "topCorner:" << topCorner.x << " " << topCorner.y << std::endl;
	std::cout << "botCorner:" << botCorner.x << " " << botCorner.y << std::endl;
	std::cout << "rightCorner:" << rightCorner.x << " " << rightCorner.y << std::endl << std::endl;


	//Replacing the missing piece with the transformed piece
	for (int i = y_min; i < y_max + 1; i++) {
	for (int j = x_min; j < x_max + 1; j++) {
	for (int k = 0; k < bytesPerPixel; k++) {
	int i_ind = 43 + i - y_min;
	int j_ind = 57 + j - x_min;
	imageData[(i * imageBreadth + j) * bytesPerPixel + k] = fitData[(i_ind * pieceBreadth + j_ind ) * bytesPerPixel +  k];
	}
	}
	}

	for (int i = y_min; i < y_max + 1; i++) {
	for (int j = x_min; j < x_max + 1; j++) {
	for (int k = 0; k < bytesPerPixel; k++) {
	int i_ind = 43 + i - y_min;
	int j_ind = 57 + j - x_min;
	outImageData[(i * imageBreadth + j) * bytesPerPixel + k] = imageData[(i * imageBreadth + j) * bytesPerPixel + k];
	}
	}
	}
	writeImage(&outImageData[0], 512, 512);*/
}

void myImage::project()
{
	int pHeight = 184; int pWidth = 350;

	// Calculated offline
	
	double H[3][3] = { 2.14872571830529,
		1.90403435442398,
		325,
		0.636802121097343,
		0.330673190459212,
		619,
		0.00278927392611717,
		0.000742710210338199, 1
	};


	vector<unsigned char>pieceImageBinary = std::vector<unsigned char>(pieceHeight * pieceBreadth);

	for (row_it = 0; row_it < pieceHeight; row_it++)
	{
		for (col_it = 0; col_it < pieceBreadth; col_it++)
		{

			pieceImageBinary[row_it * pieceBreadth + col_it] = (pieceImage[(row_it * pieceBreadth + col_it) * 3] + pieceImage[(row_it * pieceBreadth + col_it) * 3 + 1] + pieceImage[(row_it * pieceBreadth + col_it) * 3 + 2]) / 3;
		}
	}


	double xIndNew = 0, yIndNew = 0, xIndNext, yIndNext, xFloor, yFloor;
	int rowNext, colNext;
	double xCart, yCart, xCartNew, yCartNew = 0;
	unsigned char value, pixel_value;
	int debugCount = 0;
	count = 0;
	for (row_it = 0; row_it < 146; row_it++)
	{
		for (col_it = 0; col_it < 350; col_it++)
		{
			pixel_value = pieceImageBinary[(row_it * pieceBreadth + col_it)];
			//pixel_value = (pieceImage[(row_it * imageBreadth + col_it) * bytesPerPixel + 0] + pieceImage[(row_it * imageBreadth + col_it) * bytesPerPixel + 1] + pieceImage[(row_it * imageBreadth + col_it) * bytesPerPixel + 2]) / 3;
			if (pixel_value < unsigned char(250))
			{
				xIndNew = (H[0][0] * col_it + H[0][1] * row_it + H[0][2]) / (H[2][0] * col_it + H[2][1] * row_it + 1);
				yIndNew = (H[1][0] * col_it + H[1][1] * row_it + H[1][2]) / (H[2][0] * col_it + H[2][1] * row_it + 1);
				xFloor = floor(xIndNew);
				yFloor = floor(yIndNew);
				if (xFloor > 0 && xFloor < imageBreadth - 1 && yFloor > 0 && yFloor < imageHeight - 1)
				{
					for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
					{
						//count++;
						//out << count << endl;
						value = pieceImage[(row_it * pieceBreadth + col_it) * bytesPerPixel + pixel_it];
						//out << "Indx " << round(xIndNew) << " " << round(yIndNew) << "  " << col_it << "  " << row_it << "  " << value << endl;
						imageData[(yFloor * imageBreadth + xFloor) * bytesPerPixel + pixel_it] = value;
					}
				}
			}
		}
	}

	writeImage(&imageData[0], "out_trojans.raw", 972, 648, 3);


}


void myImage::ditheringMatix()
{
	// Normalize data so the matrix can be applied
	nomalize();
	ofstream out("ditherdebug.txt");
	double trial; 
	double norm;
	int N = 4;
	double fromWhere; 
	unsigned char chartrial;
	out << "normalizedData[valFrom] " << " Threshold " << "Unsigned char o/p" << "Outputdata " << endl;
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
				valFrom = (row_it * imageBreadth + col_it) ; // grey scale dithering
				valTo = valFrom;
				outImageData[valTo] =  normalizedData[valFrom] > (((double)ditherA4[((row_it * N) % N) + (col_it % N)] + 0.5) / pow(N, 2)) ? 255 : 0;
				//[REMOVE] : For debug
				/*norm = normalizedData[valFrom];
				trial = (((double)dither4[((row_it * N) % N) + (col_it % N)] + 0.5) / 16);
				fromWhere = (((row_it * N) % N) + (col_it % N)+ 0.5) / 16;
				chartrial = normalizedData[valFrom] >(((double)dither4[((row_it * N) % N) + (col_it % N)] + 0.5) / pow(N, 2)) ? 255 : 0;
				out << norm << "  " << trial << " "  << fromWhere << " " << chartrial << " " << outImageData[valTo] << endl;*/
		}
	}
}

void myImage::quantizeFourLevels()
{
	ofstream out("outQauntize.txt");
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			valFrom = row_it * imageBreadth + col_it;
			valTo = valFrom;
			outImageData[valTo] = imageData[valFrom] >= 213 ? 255 : imageData[valFrom] >= 128 && imageData[valFrom] < 213 ? (unsigned char)170 : imageData[valFrom] >= 43 ? (unsigned char)85 : 0;
			/* [REMOVE] For Debug
			out << valTo << " " << int(imageData[valFrom]) << "  " << int(outImageData[valTo]) << endl;*/
		}
	}
}

void myImage::errorDiffuseFS()
{
	ofstream out("errorFS.txt");
	float error = 0;
	int count = 0;
	nomalize();
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		if (row_it % 2 == 0)
		{
			for (col_it = 0; col_it < imageBreadth; col_it++)
			{
				valFrom = row_it * imageBreadth + col_it;
				if (normalizedData[valFrom] >= 0.5)
				{
					error = normalizedData[valFrom] - 1;
					normalizedData[valFrom] = 1;
				}
				else
				{
					error = normalizedData[valFrom] - 0;
					normalizedData[valFrom] = 0;
				}
				 
				if (col_it < imageBreadth -1 && row_it < imageHeight  -1)
				{
					count++;
					if((row_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + col_it] = (double)5 / 16 * error + normalizedData[(row_it + 1) * imageBreadth + col_it];
					if((row_it + 1) < (imageHeight - 1) && col_it - 1 > 0) normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)] = (double)3 / 16 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)];
					if((row_it + 1) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)] = (double)1 / 16 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)];
					if((col_it + 1) < (imageBreadth - 1)) normalizedData[row_it * imageBreadth + (col_it + 1)] =  (double)7 / 16 * error + normalizedData[row_it * imageBreadth + (col_it + 1)];
					out << count << endl;
				}
		}
			
		}


		else
		{
			for (col_it = imageBreadth - 1 ; col_it >= 0; col_it--)
			{
				valFrom = row_it * imageBreadth + col_it;
				if (normalizedData[valFrom] >= 0.5)
				{
					error = normalizedData[valFrom] - 1;
					normalizedData[valFrom] = 1;
				}
				else
				{
					error = normalizedData[valFrom] - 0;
					normalizedData[valFrom] = 0;
				}
				if (row_it < imageHeight - 1 && col_it > 1)
				{
					if ((row_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + col_it] = (double)5 / 16 * error + normalizedData[(row_it + 1) * imageBreadth + col_it];
					if ((row_it + 1) < (imageHeight - 1) && col_it - 1 > 0) normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)] = (double)3 / 16 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)];
					if ((row_it + 1) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)] = (double)1 / 16 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)];
					if ((col_it + 1) < (imageBreadth - 1)) normalizedData[row_it * imageBreadth + (col_it + 1)] = (double)7 / 16 * error + normalizedData[row_it * imageBreadth + (col_it + 1)];

				}

			}
		}

	}

	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			valFrom = row_it * imageBreadth + col_it;
			valTo = valFrom;
			//normalizedData[valFrom] = 1 - normalizedData[valFrom];
			outImageData[valTo] = (unsigned char)round(normalizedData[valFrom] * 255);
			out << valFrom << ""  << outImageData[valTo] << endl;
		}
	}
}

void myImage::errorDiffusionJJN()
{
	ofstream out("errorFS.txt");
	float error = 0;
	int count = 0;
	nomalize();
	double rt = 0;
	rt -= 0;
	
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		if (row_it % 2 == 0)
		{
			for (col_it = 0; col_it < imageBreadth; col_it++)
			{
				
				valFrom = row_it * imageBreadth + col_it;
				if (normalizedData[valFrom] >= 0.5)
				{
					error = normalizedData[valFrom] - 1;
					normalizedData[valFrom] = 1;
				}
				else
				{
					error = normalizedData[valFrom] - 0;
					normalizedData[valFrom] = 0;
				}
				
				if (col_it < imageBreadth - 2 && row_it < imageHeight - 2)
				{
					count++;
					if ((row_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + col_it] = (double)7 / 48 * error + normalizedData[(row_it + 1) * imageBreadth + col_it];
					if ((row_it + 1) < (imageHeight - 1) && col_it - 1 > 0) normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)] = (double)5 / 48 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)] ;
					if ((row_it + 1) < (imageHeight - 1) && col_it - 2 > 0) normalizedData[(row_it + 1) * imageBreadth + (col_it - 2)] = (double)3 / 48 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it - 2)];
					if ((row_it + 1) < (imageHeight - 1) && (col_it + 2) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + (col_it + 2)] = (double)3 / 48 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it + 2)];
					if ((row_it + 1) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)] = (double)5 / 48 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)];
					
					
					if ((col_it + 1) < (imageBreadth - 1)) normalizedData[row_it * imageBreadth + (col_it + 1)] = (double)7 / 48 * error + normalizedData[row_it * imageBreadth + (col_it + 1)];
					if ((col_it + 2) < (imageBreadth - 1)) normalizedData[row_it * imageBreadth + (col_it + 2)] = (double)5 / 48 * error + normalizedData[row_it * imageBreadth + (col_it + 2)];
					
					if ((row_it + 2) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + col_it] = (double)5 / 48 * error + normalizedData[(row_it + 2) * imageBreadth + col_it];
					if ((row_it + 2) < (imageHeight - 1) && col_it - 1 > 0) normalizedData[(row_it + 2) * imageBreadth + (col_it - 1)] = (double)3 / 48 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it - 1)];
					if ((row_it + 2) < (imageHeight - 1) && col_it - 2 > 0) normalizedData[(row_it + 2) * imageBreadth + (col_it - 2)] = (double)1 / 48 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it - 2)];
					if ((row_it + 2) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + (col_it + 2)] = (double)1 / 48 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it + 2)];
					if ((row_it + 2) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + (col_it + 1)] = (double)3 / 48 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it + 1)];
				}
			}

		}
		else
		{
			for (col_it = imageBreadth - 1; col_it >= 0; col_it--)
			{
				valFrom = row_it * imageBreadth + col_it;
				if (normalizedData[valFrom] >= 0.5)
				{
					error = normalizedData[valFrom] - 1;
					normalizedData[valFrom] = 1;
				}
				else
				{
					error = normalizedData[valFrom] - 0;
					normalizedData[valFrom] = 0;
				}
				if (row_it < imageHeight - 2 && col_it > 2)
				{
					if ((row_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + col_it] = (double)7 / 48 * error + normalizedData[(row_it + 1) * imageBreadth + col_it];
					if ((row_it + 1) < (imageHeight - 1) && col_it - 1 > 0) normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)] = (double)5 / 48 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)];
					if ((row_it + 1) < (imageHeight - 1) && col_it - 2 > 0) normalizedData[(row_it + 1) * imageBreadth + (col_it - 2)] = (double)3 / 48 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it - 2)];
					if ((row_it + 1) < (imageHeight - 1) && (col_it + 2) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + (col_it + 2)] = (double)3 / 48 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it + 2)];
					if ((row_it + 1) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)] = (double)5 / 48 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)];


					if ((col_it + 1) < (imageBreadth - 1)) normalizedData[row_it * imageBreadth + (col_it + 1)] = (double)7 / 48 * error + normalizedData[row_it * imageBreadth + (col_it + 1)];
					if ((col_it + 2) < (imageBreadth - 1)) normalizedData[row_it * imageBreadth + (col_it + 2)] = (double)5 / 48 * error + normalizedData[row_it * imageBreadth + (col_it + 2)];

					if ((row_it + 2) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + col_it] = (double)5 / 48 * error + normalizedData[(row_it + 2) * imageBreadth + col_it];
					if ((row_it + 2) < (imageHeight - 1) && col_it - 1 > 0) normalizedData[(row_it + 2) * imageBreadth + (col_it - 1)] = (double)3 / 48 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it - 1)];
					if ((row_it + 2) < (imageHeight - 1) && col_it - 2 > 0) normalizedData[(row_it + 2) * imageBreadth + (col_it - 2)] = (double)1 / 48 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it - 2)];
					if ((row_it + 2) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + (col_it + 2)] = (double)1 / 48 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it + 2)];
					if ((row_it + 2) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + (col_it + 1)] = (double)3 / 48 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it + 1)];
				}

			}
		}

	}

	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{

			
			valFrom = row_it * imageBreadth + col_it;
			valTo = valFrom;
			//normalizedData[valFrom] = 1 - normalizedData[valFrom];
			outImageData[valTo] = (unsigned char)round(normalizedData[valFrom] * 255);
			out << valFrom << "" << outImageData[valTo] << endl;
		}
	}
}

void myImage::errorDiffusionStucki()
{
	ofstream out("errorFS.txt");
	float error = 0;
	int count = 0;
	nomalize();
	double rt = 0;
	rt -= 0;

	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		if (row_it % 2 == 0)
		{
			for (col_it = 0; col_it < imageBreadth; col_it++)
			{

				valFrom = row_it * imageBreadth + col_it;
				if (normalizedData[valFrom] >= 0.5)
				{
					error = normalizedData[valFrom] - 1;
					normalizedData[valFrom] = 1;
				}
				else
				{
					error = normalizedData[valFrom] - 0;
					normalizedData[valFrom] = 0;
				}

				if (col_it < imageBreadth - 2 && row_it < imageHeight - 2)
				{
					count++;
					if ((row_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + col_it] = (double)8 /42 * error + normalizedData[(row_it + 1) * imageBreadth + col_it];
					if ((row_it + 1) < (imageHeight - 1) && col_it - 1 > 0) normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)] = (double)4/ 42 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)];
					if ((row_it + 1) < (imageHeight - 1) && col_it - 2 > 0) normalizedData[(row_it + 1) * imageBreadth + (col_it - 2)] = (double)2/ 42 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it - 2)];
					if ((row_it + 1) < (imageHeight - 1) && (col_it + 2) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + (col_it + 2)] = (double)2/ 42 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it + 2)];
					if ((row_it + 1) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)] = (double)4/ 42 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)];


					if ((col_it + 1) < (imageBreadth - 1)) normalizedData[row_it * imageBreadth + (col_it + 1)] = (double)8 / 42 * error + normalizedData[row_it * imageBreadth + (col_it + 1)];
					if ((col_it + 2) < (imageBreadth - 1)) normalizedData[row_it * imageBreadth + (col_it + 2)] = (double)4 / 42 * error + normalizedData[row_it * imageBreadth + (col_it + 2)];

					if ((row_it + 2) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + col_it] = (double)4/ 42* error + normalizedData[(row_it + 2) * imageBreadth + col_it];
					if ((row_it + 2) < (imageHeight - 1) && col_it - 1 > 0) normalizedData[(row_it + 2) * imageBreadth + (col_it - 1)] = (double)2/ 42 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it - 1)];
					if ((row_it + 2) < (imageHeight - 1) && col_it - 2 > 0) normalizedData[(row_it + 2) * imageBreadth + (col_it - 2)] = (double)1 / 42 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it - 2)];
					if ((row_it + 2) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + (col_it + 2)] = (double)1 / 42 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it + 2)];
					if ((row_it + 2) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + (col_it + 1)] = (double)2 / 42 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it + 1)];
				}
			}

		}
		else
		{
			for (col_it = imageBreadth - 1; col_it >= 0; col_it--)
			{
				valFrom = row_it * imageBreadth + col_it;
				if (normalizedData[valFrom] >= 0.5)
				{
					error = normalizedData[valFrom] - 1;
					normalizedData[valFrom] = 1;
				}
				else
				{
					error = normalizedData[valFrom] - 0;
					normalizedData[valFrom] = 0;
				}
				if (row_it < imageHeight - 2 && col_it > 2)
				{
					if ((row_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + col_it] = (double)8 / 42 * error + normalizedData[(row_it + 1) * imageBreadth + col_it];
					if ((row_it + 1) < (imageHeight - 1) && col_it - 1 > 0) normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)] = (double)4 / 42 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it - 1)];
					if ((row_it + 1) < (imageHeight - 1) && col_it - 2 > 0) normalizedData[(row_it + 1) * imageBreadth + (col_it - 2)] = (double)2 / 42 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it - 2)];
					if ((row_it + 1) < (imageHeight - 1) && (col_it + 2) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + (col_it + 2)] = (double)2 / 42 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it + 2)];
					if ((row_it + 1) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)] = (double)4 / 42 * error + normalizedData[(row_it + 1) * imageBreadth + (col_it + 1)];


					if ((col_it + 1) < (imageBreadth - 1)) normalizedData[row_it * imageBreadth + (col_it + 1)] = (double)8 / 42 * error + normalizedData[row_it * imageBreadth + (col_it + 1)];
					if ((col_it + 2) < (imageBreadth - 1)) normalizedData[row_it * imageBreadth + (col_it + 2)] = (double)4 / 42 * error + normalizedData[row_it * imageBreadth + (col_it + 2)];

					if ((row_it + 2) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + col_it] = (double)4 / 42 * error + normalizedData[(row_it + 2) * imageBreadth + col_it];
					if ((row_it + 2) < (imageHeight - 1) && col_it - 1 > 0) normalizedData[(row_it + 2) * imageBreadth + (col_it - 1)] = (double)2 / 42 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it - 1)];
					if ((row_it + 2) < (imageHeight - 1) && col_it - 2 > 0) normalizedData[(row_it + 2) * imageBreadth + (col_it - 2)] = (double)1 / 42 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it - 2)];
					if ((row_it + 2) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + (col_it + 2)] = (double)1 / 42 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it + 2)];
					if ((row_it + 2) < (imageHeight - 1) && (col_it + 1) < (imageBreadth - 1)) normalizedData[(row_it + 2) * imageBreadth + (col_it + 1)] = (double)2 / 42 * error + normalizedData[(row_it + 2) * imageBreadth + (col_it + 1)];
				}

			}
		}

	}

	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{


			valFrom = row_it * imageBreadth + col_it;
			valTo = valFrom;
			//normalizedData[valFrom] = 1 - normalizedData[valFrom];
			outImageData[valTo] = (unsigned char)round(normalizedData[valFrom] * 255);
			out << valFrom << "" << outImageData[valTo] << endl;
		}
	}
}

void myImage::convertToBinary(unsigned char imageData[], unsigned char imageDataBinary[])
{
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			valFrom = row_it * imageBreadth + col_it;
			valFromR = valFrom * bytesPerPixel + 0;
			valFromG = valFrom * bytesPerPixel + 1;
			valFromB = valFrom * bytesPerPixel + 2;
			valTo = valFrom;
			imageDataBinary[valTo] = (unsigned char)((imageData[valFromR] + imageData[valFromR] + imageData[valFromR]) / 3);
		}
	}
}

void myImage::binarize(unsigned char imageDataBinary[], unsigned char imageDataOut[])
{
	double sum = 0;
	unsigned char cmp;
	int count = 0;
	int curRow, curCol;
	int check = 0;
	for (int y = 0; y < imageHeight; ++y) {
		for (int x = 0; x < imageBreadth; ++x)
		{
			valFrom = y * imageBreadth + x;
			valTo = valFrom;
			sum = 0;
			for (int row_fil = 0; row_fil < 35; row_fil++)
			{
			for (int col_fil = 0; col_fil < 35; col_fil++)
			{
			curRow = y + row_fil - 17;
			curCol = x + col_fil - 17;
			if (curRow >= 0 && curRow <= imageHeight -1 && curCol >= 0 && curCol <= imageBreadth - 1)
			{

			count++;
			//out << curRow << " " << curCol <<  "  " << count << " " << check << endl;
			sum += imageDataBinary[(curRow)* imageBreadth + (curCol)];
			}
			}
			}
			check++;
			//out << check << " " << count << endl;
			cmp = (unsigned char)(sum / count);
			imageDataOut[valTo] = (imageDataBinary[valFrom] > cmp) ? 255 : 0;
			count = 0; 
		}
	}

	ofstream outImg("img_25.txt");
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			outImg << imageDataOut[row_it * imageBreadth + col_it];
		}
	}
	
	/*	ifstream in("img_91.txt");
	for (int i = 0; i < imageBreadth * imageHeight; i++)
	{
		in >> imageDataOut[i];

	}
	*/

}


void findMode(unsigned char arr[], int len)
{
	int count = 0;
	for (int i = 0; i < len; i++)
	{
		if (arr[i] == 0)count++;
	}

	if (count >= len / 2)
	{
		for (int i = 0; i < len; i++)
		{
			if (arr[i] = 0);
		}
	}
	else
	{
		for (int i = 0; i < len; i++)
		{
			if (arr[i] = 1);
		}
	}

}



void myImage::medianFilter(unsigned char imageDataBinary[], unsigned char imageDataOut[], int N)
{
	int curRow = 0, curCol = 0;
	int size = 2 * N + 1;
	vector<unsigned char> findMedianOf(size * size);
	vector<unsigned char> deNoiseImage(imageBreadth * imageHeight);
	unsigned char med;
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			count = 0;
			for (int row_fil = 0; row_fil < size; row_fil++)
			{
				for (int col_fil = 0; col_fil < size; col_fil++)
				{
					curRow = row_it + row_fil - N;
					curCol = col_it + col_fil - N;
					if (curRow > 0 && curRow < imageHeight - 1 && curCol > 0 && curCol < imageBreadth - 1)
					{
						findMedianOf[count++] = imageDataBinary[curRow * imageBreadth + curCol];
					}
				}
			}
			sort(findMedianOf.begin(), findMedianOf.end());
			med = findMedianOf[count / 2];
			imageDataOut[row_it * imageBreadth + col_it] = (unsigned char)med;
		}
	}
	/*
	vector<unsigned char> mode(5*5);
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			count = 0;
			for (int row_fil = 0; row_fil < 3; row_fil++)
			{
				for (int col_fil = 0; col_fil < 3; col_fil++)
				{
					curRow = row_it + row_fil - 1;
					curCol = col_it + col_fil - 1;
					if (curRow > 0 && curRow < imageHeight - 1 && curCol > 0 && curCol < imageBreadth - 1)
					{
						mode[count++] = imageDataOut[curRow * imageBreadth + curCol];
					}
				}
			}
			findMode(&mode[0], count);
			int it = count;
			count = 0;
			for (int row_fil = 0; row_fil < 5; row_fil++)
			{
				for (int col_fil = 0; col_fil < 5; col_fil++)
				{
					curRow = row_it + row_fil - 2;
					curCol = col_it + col_fil - 2;
					if (curRow > 0 && curRow < imageHeight - 1 && curCol > 0 && curCol < imageBreadth - 1)
					{
						if (count <= it)
						{
							imageDataOut[curRow * imageBreadth + curCol] = mode[count++];
						}
						 
					}
				}
			}*/
		
	

}

void myImage::convertFromDisplay(unsigned char imageDataOut[])
{
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			imageDataOut[row_it * imageBreadth + col_it] = (imageDataOut[row_it * imageBreadth + col_it] == 255) ? 1 : 0;
		}
	}

}

int myImage::convertToDisplay(unsigned char imageDataOut[])
{
	int countGrain = 0;
	int filterSize = 25;
	ofstream out("graincount.txt");
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			valFrom = row_it * imageBreadth + col_it;
			if (imageDataOut[valFrom] == (unsigned char)(1))
			{
				imageDataOut[valFrom] = 255;
				countGrain++;

			}
		}

	}
	cout << countGrain / filterSize << endl;
	out << countGrain / filterSize << endl;
	return countGrain / filterSize;
}
