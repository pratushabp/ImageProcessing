
#include "myImage.hpp"
#include "matrix.h"
#include <iostream>
#include <math.h>

#define M_PI 3.141592653589793

ofstream rankOut("rank.txt");

struct coord {
	int x = 0;
	int y = 0;
};

void myImage::readImage( char inFile[], int imageHeight, int imageBreadth, int bytesPerPixel)
{
	FILE* imageRead = 0;
	if (fopen_s(&imageRead, inFile, "rb"))
	{
		cout << " cannot read " << inFile << endl;
	}
	fread(&imageData[0], sizeof(unsigned char), imageHeight * imageBreadth * bytesPerPixel, imageRead);




}


void myImage::writeImage(unsigned char toWrite[], char* outFile, int height, int breadth, int bytes)
{

	FILE* outImage = 0;
	if (fopen_s(&outImage, outFile, "wb"))
	{
		cout << "Cannot open input file" << endl;
		exit(0);
	}
	fwrite(toWrite, sizeof(unsigned char), height * breadth * bytes, outImage);
}

double* tensor(double transposed[], double untransposed[], double tensorProduct[])
{
	int size = 5;
	for (int i = 0; i < size; i++)
	{
		int c = 0;
		for (int j = 0; j < size; j++)
		{

			tensorProduct[i * size + j] = transposed[i] * untransposed[c++];
		}
	}
	return tensorProduct;
}



vector<unsigned char> myImage::applyFilter(double filter[])
{
	vector<unsigned char> temp = std::vector<unsigned char>((imageHeight + (2 * imageExtend)) * (imageBreadth + (2 * imageExtend)) * bytesPerPixel);
	//return temp;
	int exImageBreadth = imageBreadth + (2 * imageExtend);
	int exImageHeight = imageHeight + (2 * imageExtend);
	int pixel_count = 0;
	double sum = 0;
	vector<unsigned char>window(25);

	for (row_it = imageExtend; row_it < imageExtend + imageHeight; row_it++)
	{
		for (col_it = imageExtend; col_it < imageExtend + imageHeight; col_it++)
		{
			pixel_count = 0;
			sum = 0;
			for (int winRow = row_it - N; winRow <= row_it + N; winRow++)
			{
				for (int winCol = col_it - N; winCol <= col_it + N; winCol++)
				{
					if (winRow < exImageHeight && winRow >= 0 && winCol < exImageBreadth && winCol >= 0)
					{
						sum += filter[pixel_count] * boundaryImage[winRow * exImageBreadth + winCol];
						pixel_count++;

					}
				}
			}
			//featureVector
			temp[row_it * exImageBreadth + col_it] = sum / pixel_count;

		}
	}

	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			valFrom = (row_it + imageExtend)* exImageBreadth + (col_it + imageExtend);
			valTo = row_it  * imageBreadth + col_it;
			imageData[valTo] = temp[valFrom];
		}
	}
	Mat lawImg;
	lawImg.create(cvSize(imageBreadth, imageHeight), CV_8UC1);
	memcpy(lawImg.data, &imageData[0], 128 * 128);
	imwrite("law_1.jpg", lawImg);

	return imageData;

}

vector<double> myImage::textClassify()
{
	

	int exImageBreadth = imageBreadth + (2 * imageExtend);
	int exImageHeight = imageHeight + (2 * imageExtend);

	double filters[5][5] = { { 1.0,  4.0,  6.0,  4.0,  1.0 } ,{ -1.0, -2.0 , 0.0,  2.0,  1.0 } ,{ -1.0,  0.0,  2.0,  0.0, -1.0 } ,{ -1.0,  2.0,  0.0, -2.0,  1.0 } ,{ 1.0,  -4.0,  6.0, -4.0,  1.0 } };

	double L5[5] = { 1.0,  4.0,  6.0,  4.0,  1.0 };
	int E5[5] = { -1.0, -2.0 , 0.0,  2.0,  1.0 };
	int S5[5] = { -1.0,  0.0,  2.0,  0.0, -1.0 };
	int W5[5] = { -1.0,  2.0,  0.0, -2.0,  1.0 };
	int R5[5] = { 1.0,  -4.0,  6.0, -4.0,  1.0 };
	vector<double>product(25);
	double what[25];
	vector<double>product1(25);
	int pixel_count = 0;
	int featureCount = 0;
	int sum = 0;
	vector<vector<double>> tensorProduct;
	double lol[25];
	tensor(L5, L5, lol);
	Mat img, imgFil, kernel;
	img.create(cvSize(128, 128), CV_8UC1);
	kernel.create(cvSize(5, 5), CV_32F);
	memcpy(img.data, &imageData[0], 128 * 128);
	memcpy(kernel.data, &lol[0], 25);

	
	filter2D(img, imgFil, -1, kernel, Point(-1, -1), 0, BORDER_DEFAULT);
	imgFil = cv::abs(imgFil);
	imwrite("please.jpg", imgFil);
	
	vector<vector<unsigned char>> filteredImages;
	vector<double>featureVector;
	vector<unsigned char> temp = vector<unsigned char>(128 * 128 * 1);
	//std::vector<unsigned char>im = std::vector<unsigned char>((imageHeight + (2 * imageExtend)) * (imageBreadth + (2 * imageExtend)) * bytesPerPixel);

	
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			tensor(filters[i], filters[j], &product[0]);
			tensorProduct.push_back(product);
		}
	}
	
	for (int tensorIt = 0; tensorIt < tensorProduct.size(); tensorIt++)
	{
		product1 = (tensorProduct[tensorIt]);
		temp = applyFilter(&product1[0]);
		filteredImages.push_back(temp);
	}

	double energyOfFilter = 0;
	for (int it = 0; it < filteredImages.size(); it++)
	{
		temp = filteredImages[it];
		energyOfFilter = 0;
		for (row_it = 0; row_it < imageHeight; row_it++)
		{
			for (col_it = 0; col_it < imageBreadth; col_it++)
			{
				energyOfFilter += temp[row_it * imageBreadth + col_it];
			}
		}
			featureVector.push_back(energyOfFilter / (imageHeight * imageBreadth));
	}
	

	return featureVector;

	Mat lawImg;
	lawImg.create(cvSize(imageBreadth, imageHeight), CV_8UC1);
	memcpy(lawImg.data, &temp[0], 128 * 128);
	//imwrite("law_1.jpg", lawImg);
	
}


vector<vector<unsigned char>> myImage::textSegment()
{


	int exImageBreadth = imageBreadth + (2 * imageExtend);
	int exImageHeight = imageHeight + (2 * imageExtend);

	double filters[5][5] = { { 1.0,  4.0,  6.0,  4.0,  1.0 } ,{ -1.0, -2.0 , 0.0,  2.0,  1.0 } ,{ -1.0,  0.0,  2.0,  0.0, -1.0 } ,{ -1.0,  2.0,  0.0, -2.0,  1.0 } ,{ 1.0,  -4.0,  6.0, -4.0,  1.0 } };

	double L5[5] = { 1.0,  4.0,  6.0,  4.0,  1.0 };
	int E5[5] = { -1.0, -2.0 , 0.0,  2.0,  1.0 };
	int S5[5] = { -1.0,  0.0,  2.0,  0.0, -1.0 };
	int W5[5] = { -1.0,  2.0,  0.0, -2.0,  1.0 };
	int R5[5] = { 1.0,  -4.0,  6.0, -4.0,  1.0 };
	vector<double>product(25);
	double what[25];
	vector<double>product1(25);
	int pixel_count = 0;
	int featureCount = 0;
	int sum = 0;
	vector<vector<double>> tensorProduct;
	double lol[25];
	tensor(L5, L5, lol);
	Mat img, imgFil, kernel;
	img.create(cvSize(128, 128), CV_8UC1);
	kernel.create(cvSize(5, 5), CV_32F);
	//memcpy(img.data, &imageData[0], 128 * 128);
	memcpy(kernel.data, &lol[0], 25);


	filter2D(img, imgFil, -1, kernel, Point(-1, -1), 0, BORDER_DEFAULT);
	imgFil = cv::abs(imgFil);
	imwrite("please.jpg", imgFil);

	vector<vector<unsigned char>> filteredImages;
	vector<vector<unsigned char>>featureVector;
	vector<unsigned char> temp = vector<unsigned char>(imageBreadth * imageHeight * 1);
	vector<unsigned char> tempLocalMean = vector<unsigned char>(imageBreadth * imageHeight * 1);
	//std::vector<unsigned char>im = std::vector<unsigned char>((imageHeight + (2 * imageExtend)) * (imageBreadth + (2 * imageExtend)) * bytesPerPixel);


	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			tensor(filters[i], filters[j], &product[0]);
			tensorProduct.push_back(product);
		}
	}

	for (int tensorIt = 0; tensorIt < tensorProduct.size(); tensorIt++)
	{
		product1 = (tensorProduct[tensorIt]);
		temp = applyFilter(&product1[0]);
		filteredImages.push_back(temp);
	}

	double energyOfFilter = 0;
	for (int it = 0; it < filteredImages.size(); it++)
	{
		temp = filteredImages[it];
		energyOfFilter = 0;
		for (row_it = 0; row_it < imageHeight; row_it++)
		{
			for (col_it = 0; col_it < imageBreadth; col_it++)
			{
				sum = 0;
				pixel_count = 0;
				for (int winRow = row_it - N; winRow <= row_it + N; winRow++)
				{
					for (int winCol = col_it - N; winCol <= col_it + N; winCol++)
					{
						if (winRow < imageHeight && winRow >= 0 && winCol < imageBreadth && winCol >= 0)
						{
							//window[pixel_count] = boundaryImage[winRow * 5 + winCol];
							sum += temp[winRow * imageBreadth + winCol];
							pixel_count++;

						}
					}
				}

				tempLocalMean[row_it * imageBreadth + col_it] = (sum / pixel_count);
				//energyOfFilter += temp[row_it * imageBreadth + col_it];
			}
		}
		featureVector.push_back(tempLocalMean);
	}

	unsigned char m;// = featureVector[0][0];
	vector<unsigned char>pixelVector;
	vector<vector<unsigned char>>features;
	for (int i = 0; i < imageBreadth * imageHeight; i++)
	{
		for (int j = 0; j < featureVector.size(); j++)
		{
			m = (featureVector[j][i]);
			pixelVector.push_back(m);
		}
		features.push_back(pixelVector);
		pixelVector.clear();
	}


	return features;

	Mat lawImg;
	lawImg.create(cvSize(imageBreadth, imageHeight), CV_8UC1);
	//memcpy(lawImg.data, &temp[0], 128 * 128);
	//imwrite("law_1.jpg", lawImg);

}

vector<vector<unsigned char>> myImage::adtextSegment()

{


	int exImageBreadth = imageBreadth + (2 * imageExtend);
	int exImageHeight = imageHeight + (2 * imageExtend);

	double filters[5][5] = { { 1.0,  4.0,  6.0,  4.0,  1.0 } ,{ -1.0, -2.0 , 0.0,  2.0,  1.0 } ,{ -1.0,  0.0,  2.0,  0.0, -1.0 } ,{ -1.0,  2.0,  0.0, -2.0,  1.0 } ,{ 1.0,  -4.0,  6.0, -4.0,  1.0 } };

	double L5[5] = { 1.0,  4.0,  6.0,  4.0,  1.0 };
	int E5[5] = { -1.0, -2.0 , 0.0,  2.0,  1.0 };
	int S5[5] = { -1.0,  0.0,  2.0,  0.0, -1.0 };
	int W5[5] = { -1.0,  2.0,  0.0, -2.0,  1.0 };
	int R5[5] = { 1.0,  -4.0,  6.0, -4.0,  1.0 };
	vector<double>product(25);
	double what[25];
	vector<double>product1(25);
	int pixel_count = 0;
	int featureCount = 0;
	int sum = 0;
	vector<vector<double>> tensorProduct;
	double lol[25];
	tensor(L5, L5, lol);
	Mat img, imgFil, kernel;
	img.create(cvSize(128, 128), CV_8UC1);
	kernel.create(cvSize(5, 5), CV_32F);
	//memcpy(img.data, &imageData[0], 128 * 128);
	memcpy(kernel.data, &lol[0], 25);


	filter2D(img, imgFil, -1, kernel, Point(-1, -1), 0, BORDER_DEFAULT);
	imgFil = cv::abs(imgFil);
	imwrite("please.jpg", imgFil);

	vector<vector<unsigned char>> filteredImages;
	vector<vector<unsigned char>>featureVector;
	vector<unsigned char> temp = vector<unsigned char>(imageBreadth * imageHeight * 1);
	vector<unsigned char> tempLocalMean = vector<unsigned char>(imageBreadth * imageHeight * 1);
	//std::vector<unsigned char>im = std::vector<unsigned char>((imageHeight + (2 * imageExtend)) * (imageBreadth + (2 * imageExtend)) * bytesPerPixel);


	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			tensor(filters[i], filters[j], &product[0]);
			tensorProduct.push_back(product);
		}
	}

	for (int tensorIt = 0; tensorIt < tensorProduct.size(); tensorIt++)
	{
		product1 = (tensorProduct[tensorIt]);
		temp = applyFilter(&product1[0]);
		filteredImages.push_back(temp);
	}

	double energyOfFilter = 0;
	for (int it = 0; it < filteredImages.size(); it++)
	{
		temp = filteredImages[it];
		energyOfFilter = 0;
		for (row_it = 0; row_it < imageHeight; row_it++)
		{
			for (col_it = 0; col_it < imageBreadth; col_it++)
			{
				sum = 0;
				pixel_count = 0;
				for (int winRow = row_it - N; winRow <= row_it + N; winRow++)
				{
					for (int winCol = col_it - N; winCol <= col_it + N; winCol++)
					{
						if (winRow < imageHeight && winRow >= 0 && winCol < imageBreadth && winCol >= 0)
						{
							//window[pixel_count] = boundaryImage[winRow * 5 + winCol];
							sum += temp[winRow * imageBreadth + winCol];
							pixel_count++;

						}
					}
				}
				int mean = sum / pixel_count;
				N = (mean - 0.5) * (mean +0.5); 
				sum = 0;
				pixel_count = 0;
				for (int winRow = row_it - N; winRow <= row_it + N; winRow++)
				{
					for (int winCol = col_it - N; winCol <= col_it + N; winCol++)
					{
						if (winRow < imageHeight && winRow >= 0 && winCol < imageBreadth && winCol >= 0)
						{
							//window[pixel_count] = boundaryImage[winRow * 5 + winCol];
							sum += temp[winRow * imageBreadth + winCol];
							pixel_count++;

						}
					}
				}

				tempLocalMean[row_it * imageBreadth + col_it] = (sum / pixel_count);
				//energyOfFilter += temp[row_it * imageBreadth + col_it];
			}
		}
		featureVector.push_back(tempLocalMean);
	}

	unsigned char m;// = featureVector[0][0];
	vector<unsigned char>pixelVector;
	vector<vector<unsigned char>>features;
	for (int i = 0; i < imageBreadth * imageHeight; i++)
	{
		for (int j = 0; j < featureVector.size(); j++)
		{
			m = (featureVector[j][i]);
			pixelVector.push_back(m);
		}
		features.push_back(pixelVector);
		pixelVector.clear();
	}


	return features;

	Mat lawImg;
	lawImg.create(cvSize(imageBreadth, imageHeight), CV_8UC1);
	//memcpy(lawImg.data, &temp[0], 128 * 128);
	//imwrite("law_1.jpg", lawImg);

}

void myImage::subtractMean()
{
	
	int sum = 0;
	int energyDC =  0;

	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			sum += imageData[row_it * imageBreadth + col_it];
		}
	}

	energyDC = sum / (imageBreadth * imageHeight * bytesPerPixel);


	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			subtractedImage[row_it * imageBreadth + col_it] = unsigned char(imageData[row_it * imageBreadth + col_it] - energyDC);
		}
	}

	
	//writeImage(&imageData[0], "out_mean.raw", 128, 128, 1);
	
}

void myImage::extendImage()
{
	
	int debug;
	int exImageBreadth = imageBreadth + (2 * imageExtend);
	int exImageHeight = imageHeight  +  (2 * imageExtend);

	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageBreadth; col_it++)
		{
			valFrom = row_it * imageBreadth + col_it;
			valTo = (row_it + imageExtend) * exImageBreadth + (col_it + imageExtend);
			boundaryImage[valTo] = subtractedImage[valFrom]; //(unsigned char) 255; For DebugginG
			
		}
	}
	for (row_it = 0; row_it <= imageExtend; row_it++)
	{
	for (col_it = 0; col_it < imageBreadth; col_it++)
	{	
		boundaryImage[row_it * (imageBreadth + 4) + col_it + imageExtend] = subtractedImage[imageExtend * imageBreadth + (col_it)];
		boundaryImage[(exImageHeight - 1 - row_it) * (imageBreadth + 4) + col_it + imageExtend] =  subtractedImage[(imageHeight - 1 - imageExtend) * imageBreadth + (col_it)];
		}
	}
	
	for (row_it = 0; row_it < exImageHeight ; row_it++)
	{
		imageExtend = 2;
		for (col_it = 0; col_it <= imageExtend; col_it++)
		{
			boundaryImage[row_it * exImageBreadth + col_it] = boundaryImage[row_it * (imageBreadth + 4) + col_it + imageExtend];
			boundaryImage[(row_it) * exImageBreadth + ( exImageBreadth - 1 - col_it )] = boundaryImage[(exImageHeight - 1 - row_it) * (imageBreadth + 4) + col_it + imageExtend];
		}
	}	

	
}



void myImage::cannyEdge()
{
	Mat img, imgGray, imgBlur, dummyForThresh,edgeMap;
	double threshMax, threshMin;
	Scalar mu, sigma;
	img.create(cvSize(imageBreadth, imageHeight), CV_8UC3);
	memcpy(img.data, &imageData[0], imageHeight * imageBreadth * 3);
	cvtColor(img, imgGray, CV_RGB2GRAY);
	meanStdDev(imgGray, mu, sigma);
	threshMax = mu.val[0] + 0.5 * sigma.val[0];
	threshMin = mu.val[0] - 0.5  * sigma.val[0];
	//threshMax = threshold(imgGray, dummyForThresh, 0, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
	//threshMin = 0.5 * threshMax;
	Canny(imgGray, edgeMap, threshMin, threshMax, 3, true);
	imwrite("Jaguar_mean.jpg", edgeMap);
	namedWindow("display", CV_WINDOW_AUTOSIZE);
	imshow("display", edgeMap);
	waitKey(0);
	destroyAllWindows();


	
}