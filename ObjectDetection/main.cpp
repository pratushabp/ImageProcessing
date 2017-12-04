/* EE569 Homework Assignment #3
* Date: 10 / 16 / 16
* Name : Pratusha Bhuvana Prasad
* ID : 6169305935
* Email : pratushp@usc.edu
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
#include <malloc.h>
#include "textureAnalysis.h"
#include "featureExtractionAndMatching.h"
#include "edgeDetection.h"
#include "opencv.hpp"

#define level 256
#define M_PI 3.141592653589793
using namespace::std;




int main(int argc, char* argv[])
{
	cout << "Digital Image Processing" << endl;
	cout << "Assignment 3 " << endl;
	cout << "Pratusha Bhuvana Prasad" << endl;


	if (argc < 2)
	{
		{
			cerr << "Usage:" << endl;
			cerr << "puzzleMatch" << endl;
			cerr << "homographicTransform" << endl;
			cerr << "ditheringMatrix" << endl;
			cerr << "errorDiffusion" << endl;
			cerr << "riceGrain" << endl;
			cerr << "match" << endl;
		}
	}
	else if (!(strcmp(argv[1], "texClassify")))
		texClassify(argc - 1, argv + 1);
	else if (!(strcmp(argv[1], "texSegmentation")))
		texSegmentation(argc - 1, argv + 1);
	else if (!(strcmp(argv[1], "advancedTexSegmentation")))
		advancedTexSegmentation(argc - 1, argv + 1);
	else if (!(strcmp(argv[1], "featureExtraction")))
		featureExtraction(argc - 1, argv + 1);
	else if (!(strcmp(argv[1], "imageMatch")))
		imageMatch(argc - 1, argv + 1);
	else if (!(strcmp(argv[1], "bagOfWords")))
		bagOfWords(argc - 1, argv + 1);
	else if (!(strcmp(argv[1], "cannyEdge")))
		cannyEdge(argc - 1, argv + 1);
	else if (!(strcmp(argv[1], "StructuredEdge")))
		structuredEdge(argc - 1, argv + 1);
	else if (!(strcmp(argv[1], "performanceEval")))
		performanceEval(argc - 1, argv + 1);

}