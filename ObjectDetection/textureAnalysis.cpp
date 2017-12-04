#include "textureAnalysis.h"


#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>


void texClassify(int argc, char * argv[])
{
	cout << "texture classification" << endl;
	if (argc < 5)
	{
		cout << "Usage: " << argv[0] << "  inputimage  Height Breadth bytesPerPixel" << endl;
		exit(0);
	}
	//Instantiate class
	myImage inputImage(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), argv[1]);
	vector<vector<double>>featureVector;
	vector<double>temp;


	char* files[]= { "Texture1.raw", "Texture2.raw", "Texture3.raw", "Texture4.raw", "Texture5.raw", "Texture6.raw", "Texture7.raw", "Texture8.raw","Texture9.raw", "Texture10.raw", "Texture11.raw", "Texture12.raw" };
	for (int img = 0; img < 12; img++)
	{
		inputImage.readImage(files[img], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
		inputImage.subtractMean();
		inputImage.extendImage();
		temp = inputImage.textClassify();
		featureVector.push_back(temp);
	}
	float i;
	Mat features;
	features.create(cvSize(25, 12), CV_32F);
	
	for (int row_it = 0; row_it < 12; row_it++)
	{
		for (int col_it = 0; col_it < 25; col_it++)
		{
			features.at<float>(row_it, col_it) = featureVector[row_it][col_it];
		}
	}

		Mat calcMean, reduced, labels, clusterCenter;
		PCA pca(features, calcMean, CV_PCA_DATA_AS_ROW, 25);
		pca.project(features, reduced);
		TermCriteria forKmeans = TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 300, 0.3);

		kmeans(reduced, 4, labels, forKmeans, 10, KMEANS_PP_CENTERS, clusterCenter);

		
		cout << labels << endl;
		

}


void texSegmentation(int argc, char * argv[])
{

	cout << "texture classification" << endl;
	if (argc < 5)
	{
		cout << "Usage: " << argv[0] << "  inputimage  Height Breadth bytesPerPixel" << endl;
		exit(0);
	}
	//Instantiate class
	myImage inputImage(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), argv[1]);
	vector<vector<double>>featureVector;
	vector<vector<unsigned char>>temp;


	//char* files[] = { "Texture1.raw", "Texture2.raw", "Texture3.raw", "Texture4.raw", "Texture5.raw", "Texture6.raw", "Texture7.raw", "Texture8.raw","Texture9.raw", "Texture10.raw", "Texture11.raw", "Texture12.raw" };
	//for (int img = 0; img < 12; img++)
	//{
		inputImage.readImage(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
		inputImage.subtractMean();
		inputImage.extendImage();
		temp = inputImage.textSegment();
		
		
	//}
	int i = 0;
	ofstream out("GOD.txt");
	Mat features;
	features.create(cvSize(25, inputImage.imageBreadth * inputImage.imageHeight), CV_32F);

	for (int row_it = 0; row_it < inputImage.imageBreadth * inputImage.imageHeight; row_it++)
	{
		for (int col_it = 0; col_it < 25; col_it++)
		{
			features.at<float>(row_it, col_it) = temp[row_it][col_it];
			out << i++ << endl;	
		}
	}

	Mat calcMean, reduced, labels, clusterCenter;
	PCA pca(features, calcMean, CV_PCA_DATA_AS_ROW, 1);
	pca.project(features, reduced);
	TermCriteria forKmeans = TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 300, 0.3);

	kmeans(reduced, 5, labels, forKmeans, 10, KMEANS_PP_CENTERS, clusterCenter);

	cout << labels << endl;

	vector<int>displayLabel(inputImage.imageBreadth * inputImage.imageHeight);
	
	if (labels.isContinuous())
	{
		displayLabel.assign((int*)labels.datastart, (int*)labels.dataend);
	}
	else
	{
		for (int i = 0; i < labels.rows; i++)
		{
			displayLabel.insert(displayLabel.end(), (int*)labels.ptr<unsigned char>(i), (int*)labels.ptr<unsigned char>(i) + labels.cols);

		}
	}
	vector<unsigned char>displayImg(inputImage.imageBreadth * inputImage.imageHeight);
	//	memcpy(&displayLabel[0], labels.data, inputImage.imageBreadth *inputImage.imageHeight);
	int replace;
	unsigned char value;

	for (int row_it = 0; row_it < inputImage.imageHeight; row_it++)
	{
		for (int col_it = 0; col_it < inputImage.imageBreadth; col_it++)
		{

			replace = displayLabel[row_it * inputImage.imageBreadth + col_it];
				if (replace == 0) value = 0;
				else if (replace == 1) value = 63;
				else if (replace == 2) value = 127;
				else if (replace == 3) value = 199;
				else if (replace == 4) value = 255;
				displayImg[row_it * inputImage.imageBreadth + col_it] = value;
		}
	}
	inputImage.writeImage(&displayImg[0], "segment.raw", inputImage.imageHeight, inputImage.imageBreadth, 1);
	Mat display;
	display.create(cvSize(inputImage.imageBreadth, inputImage.imageHeight), CV_8UC1);
	memcpy(display.data, &displayImg[0], inputImage.imageBreadth *inputImage.imageHeight);
	imwrite("Segmentation.jpg", display);

}

void advancedTexSegmentation(int argc, char * argv[])

{

	cout << "texture classification" << endl;
	if (argc < 5)
	{
		cout << "Usage: " << argv[0] << "  inputimage  Height Breadth bytesPerPixel" << endl;
		exit(0);
	}
	//Instantiate class
	myImage inputImage(atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), argv[1]);
	vector<vector<double>>featureVector;
	vector<vector<unsigned char>>temp;
	Mat img1;
	img1 = imread("segment.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	imwrite("segment1.jpg", img1);

	//char* files[] = { "Texture1.raw", "Texture2.raw", "Texture3.raw", "Texture4.raw", "Texture5.raw", "Texture6.raw", "Texture7.raw", "Texture8.raw","Texture9.raw", "Texture10.raw", "Texture11.raw", "Texture12.raw" };
	//for (int img = 0; img < 12; img++)
	//{
	inputImage.readImage(argv[1], atoi(argv[2]), atoi(argv[3]), atoi(argv[4]));
	inputImage.subtractMean();
	inputImage.extendImage();
	temp = inputImage.adtextSegment();


	//}
	int i = 0;
	ofstream out("GOD.txt");
	Mat features;
	features.create(cvSize(25, inputImage.imageBreadth * inputImage.imageHeight), CV_32F);

	for (int row_it = 0; row_it < inputImage.imageBreadth * inputImage.imageHeight; row_it++)
	{
		for (int col_it = 0; col_it < 25; col_it++)
		{
			features.at<float>(row_it, col_it) = temp[row_it][col_it];
			out << i++ << endl;
		}
	}

	Mat calcMean, reduced, labels, clusterCenter;
	PCA pca(features, calcMean, CV_PCA_DATA_AS_ROW, 1);
	pca.project(features, reduced);
	TermCriteria forKmeans = TermCriteria(CV_TERMCRIT_EPS + CV_TERMCRIT_ITER, 300, 0.3);

	kmeans(reduced, 5, labels, forKmeans, 10, KMEANS_PP_CENTERS, clusterCenter);

	cout << labels << endl;

	vector<int>displayLabel(inputImage.imageBreadth * inputImage.imageHeight);

	if (labels.isContinuous())
	{
		displayLabel.assign((int*)labels.datastart, (int*)labels.dataend);
	}
	else
	{
		for (int i = 0; i < labels.rows; i++)
		{
			displayLabel.insert(displayLabel.end(), (int*)labels.ptr<unsigned char>(i), (int*)labels.ptr<unsigned char>(i) + labels.cols);

		}
	}
	vector<unsigned char>displayImg(inputImage.imageBreadth * inputImage.imageHeight);
	//	memcpy(&displayLabel[0], labels.data, inputImage.imageBreadth *inputImage.imageHeight);
	int replace;
	unsigned char value;

	for (int row_it = 0; row_it < inputImage.imageHeight; row_it++)
	{
		for (int col_it = 0; col_it < inputImage.imageBreadth; col_it++)
		{

			replace = displayLabel[row_it * inputImage.imageBreadth + col_it];
			if (replace == 0) value = 0;
			else if (replace == 1) value = 63;
			else if (replace == 2) value = 127;
			else if (replace == 3) value = 199;
			else if (replace == 4) value = 255;
			displayImg[row_it * inputImage.imageBreadth + col_it] = value;
		}
	}
	inputImage.writeImage(&displayImg[0], "segment.raw", inputImage.imageHeight, inputImage.imageBreadth, 1);
	Mat display;
	display.create(cvSize(inputImage.imageBreadth, inputImage.imageHeight), CV_8UC1);
	memcpy(display.data, &displayImg[0], inputImage.imageBreadth *inputImage.imageHeight);
	imwrite("Segmentation.jpg", display);

}
