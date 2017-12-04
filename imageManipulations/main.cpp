/* EE569 Homework Assignment #2
* Date: 10/2/16
* Name: Pratusha Bhuvana Prasad
* ID: 6169305935
* Email: pratushp@usc.edu
*/

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
#include "myImage.hpp"
#include <malloc.h>

#define level 256
#define M_PI 3.141592653589793
using namespace::std;

void puzzleMatch(int argc, char* argv[])
{
			cout << "Starting puzzle Match" << endl;

		if (argc < 9)
		{
			cout << "Usage: " << argv[0] << "inputimage outimage Height Breadth bytesPerPixel pieceImage pieceImageHeight pieceImageBreadth " << endl;
			exit(0);
		}
		//Instantiate class
		myImage inputImage(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), argv[1], argv[2], atoi(argv[7]), atoi(argv[8]), argv[6]);
		inputImage.readImage();
		inputImage.pieceHeight = atoi(argv[7]);
		inputImage.pieceBreadth = atoi(argv[8]);
		inputImage.readPieceImage(atoi(argv[7]), atoi(argv[8]), atoi(argv[5]), argv[6]);
		if (!(strcmp(argv[1], "Hillary.raw")) || !(strcmp(argv[1], "hillary.raw")))
		{
			inputImage.puzzleHillary();
		}
		else
		{
			inputImage.puzzleTrump();
		}


}
void homographicTransform(int argc, char* argv[])
{
	cout << "Starting Projection" << endl;

	if (argc < 9)
	{
		cout << "Usage: " << argv[0] << "inputimage outimage Height Breadth bytesPerPixel pieceImage pieceImageHeight pieceImageBreadth " << endl;
		exit(0);
	}
	//Instantiate class
	myImage inputImage(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), argv[1], argv[2], atoi(argv[7]), atoi(argv[8]), argv[6]);
	inputImage.readImage();
	inputImage.pieceHeight = atoi(argv[7]);
	inputImage.pieceBreadth = atoi(argv[8]);
	inputImage.readPieceImage(atoi(argv[7]), atoi(argv[8]), atoi(argv[5]), argv[6]);
	inputImage.project();
	//inputImage.writeImage();

	//inputImage.writePieceImage();

}

/* Problem 2a Dithering Matrix */

void ditheringMatrix(int argc, char* argv[])
{
	cout << "Starting dithering matrix" << endl;
	if (argc < 6)
	{
		cout << "Usage: " << argv[0] << "  inputimage outimage Height Breadth bytesPerPixel" << endl;
		exit(0);
	}
	//Instantiate class
	myImage inputImage(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), argv[1], argv[2]);
	inputImage.readImage();
	//[ASK TA]: Together or separate?
	inputImage.ditheringMatix();
	//inputImage.quantizeFourLevels();
	inputImage.writeImage();

	//[FIXME] Quantize for part 2a. 3) : Close results and open again
}
void errorDiffusion(int argc, char* argv[])
{

	cout << "Starting error diffusion" << endl;
	myImage inputImage(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), argv[1], argv[2]);
	inputImage.readImage();
	//inputImage.errorDiffuseFS();
	//inputImage.errorDiffusionJJN();
	inputImage.errorDiffusionStucki();
	inputImage.writeImage();
}


int main(int argc, char* argv[])
{
	cout << "Digital Image Processing" << endl;
	cout << "Assignment 2 " << endl;
	cout << "Pratusha Bhuvana Prasad" << endl;


	if (argc < 2)
	{
		{
			cerr << "Usage:" << endl;
			cerr << "puzzleMatch" << endl;
			cerr << "homographicTransform" << endl;
			cerr << "ditheringMatrix" << endl;
			cerr << "errorDiffusion" << endl;
		}
	}

	else if (!(strcmp(argv[1], "puzzleMatch")))
		puzzleMatch(argc - 1, argv + 1);
	else if (!(strcmp(argv[1], "homographicTransform")))
		homographicTransform(argc - 1, argv + 1);
	else if (!(strcmp(argv[1], "ditheringMatrix")))
		ditheringMatrix(argc - 1, argv + 1);
	else if (!(strcmp(argv[1], "errorDiffusion")))
		errorDiffusion(argc - 1, argv + 1);

}