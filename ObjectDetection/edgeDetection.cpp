#include "edgeDetection.h"

void cannyEdge(int argc, char * argv[])
{
	Mat img;

	cout << "texture classification" << endl;
	if (argc < 5)
	{
		cout << "Usage: " << argv[0] << "  inputimage  Height Breadth bytesPerPixel" << endl;
		exit(0);
	}
	//Instantiate class
	myImage image(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), argv[1]);
	image.readImage(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
	image.cannyEdge();
	//CvMat data = cvMat(image.imageHeight, image.imageBreadth, CV_8UC3, )
	//img.create(image.imageHeight, image.imageBreadth, CV_8UC3);
	//memccpy(img.data, image.imageData, image.imageHeight * image.imageBreadth);
}

void structuredEdge(int argc, char * argv[])
{
}

void performanceEval(int argc, char * argv[])
{
}
