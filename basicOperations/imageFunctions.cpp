#include "imageFunctions.hp'"


#define level 256
#define M_PI 3.141592653589793
using namespace::std;

#define level 256
#define M_PI 3.141592653589793
using namespace::std;

struct matchPoints
{
	int pixelIntensity = 0;
	double cdf = 0.0;
};



void crop(int argc, char* argv[])
{
	FILE* imageRead = 0;
	FILE* imageWrite = 0;
	int bytesPerPixel, length, breadth;
	int rowStart, rowEnd, colStart, colEnd;
	int row_it = 0, col_it = 0, pixel_it;
	int valTo = 0, valFrom = 0;
	int cropLength = 0, cropBreadth = 0, cropSize;
	ofstream out("out.txt");

	if (argc < 9)
	{
		cerr << "Usage:" << argv[0] << "inputImage.raw" << "outputImage.raw" << "breadth" << "length" << "BytesPerPixel" << "top  left x co-ord" << "top left y co-ord" << "bottom right x co-ord" << " bottom right y co-ord" << endl;
		exit(0);
	}

	bytesPerPixel = atoi(argv[5]);
	length = atoi(argv[4]);
	breadth = atoi(argv[3]);
	rowStart = atoi(argv[6]);
	colStart = atoi(argv[7]);
	rowEnd = atoi(argv[8]);
	colEnd = atoi(argv[9]);

	if ((colEnd - colStart) < 0)
	{
		cout << "Enter parameters correctly" << endl;
		exit(0);
	}
	else
		cropLength = colEnd - colStart + 1;


	if ((rowEnd - rowStart) < 0)
	{
		cout << "Enter parameters correctly" << endl;
		exit(0);
	}
	else
		cropBreadth = rowEnd - rowStart + 1;


	vector<unsigned char> readImage(bytesPerPixel* length* breadth);
	vector<unsigned char> cropImage(bytesPerPixel* length* breadth);

	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << "Cannot open file" << argv[1] << endl;
		exit(1);
	}

	fread((unsigned char*)&readImage[0], sizeof(unsigned char), length*breadth*bytesPerPixel, imageRead);
	fclose(imageRead);
	int count = 0;
	cout << "Starting cropping" << endl;
	for (row_it = 0; row_it < cropLength - 1; ++row_it)
	{
		for (col_it = 0; col_it < cropBreadth - 1; ++col_it)
		{
			for (pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
			{
				valTo = (row_it * cropBreadth + col_it) * bytesPerPixel + pixel_it;
				valFrom = ((row_it + rowStart) * breadth + (col_it + colStart)) * bytesPerPixel + pixel_it;
				cropImage[valTo] = readImage[valFrom];
				out << count++ << endl;
			}
		}
	}
	//[FIXME] Doesn't work for (0,0) : Image iterates over if out of bounds

	if ((fopen_s(&imageWrite, argv[2], "wb")))
	{
		cout << "Cannot open output file " << endl;
	}

	cout << "Cropped image resolution(wxh) is" << cropBreadth << "x" << cropLength << endl;
	fwrite((unsigned char*)&cropImage[0], sizeof(unsigned char), cropBreadth*cropLength *bytesPerPixel, imageWrite);
	fclose(imageWrite);

}


void resize(int argc, char* argv[])
{
	int imageWidth, imageHeight, resizeWidth, resizeHeight, bytesPerPixel;
	int row_it, col_it, pixel_it;

	int valTo, valFrom;
	FILE* imageRead = 0;
	FILE* imageWrite = 0;
	ofstream out("resizeOut.txt");


	if (argc < 8)
	{
		cout << argv[0] << "input_image.raw " << "output_image.raw " << "width " << "height " << "bytesPerPixel " << "resizeWidth " << "resizeHeight " << endl;
		exit(1);
	}

	imageWidth = atoi(argv[3]);
	imageHeight = atoi(argv[4]);
	bytesPerPixel = atoi(argv[5]);
	resizeWidth = atoi(argv[6]);
	resizeHeight = atoi(argv[7]);

	vector<unsigned char> readImage(imageHeight * imageWidth * bytesPerPixel);
	vector<unsigned char> resizedImage(resizeHeight * resizeWidth * bytesPerPixel);

	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << argv[1] << endl;
		cout << "Cannot open Input image" << endl;
		exit(1);
	}
	fread(&readImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageRead);
	fclose(imageRead);

	int count = 0;
	int imageSize = imageHeight * imageWidth;
	int newSize = resizeHeight * resizeWidth;

	double resizeFactorI = double(imageHeight) / double(resizeHeight);   // Get the co-ordinates of resized image with respect to  input image
	double resizeFactorJ = double(imageWidth) / double(resizeWidth);
	for (int i = 0; i < resizeHeight - 1; i++)
	{
		for (int j = 0; j < resizeWidth - 1; j++)
		{
			for (int k = 0; k < bytesPerPixel; k++)
			{
				float i_delta = (resizeFactorI)* i;
				float j_delta = (resizeFactorJ)* j;
				int i_index = i_delta;
				int j_index = j_delta;

				int a, b, c, d;

				i_delta = i_delta - i_index;
				j_delta = j_delta - j_index;

				// forward and backward index 
				int r_index = j_index + 1;
				int b_index = i_index + 1;

				if (r_index >= imageWidth)
					r_index %= resizeWidth;


				if (b_index >= imageHeight)
					b_index %= resizeHeight;



				a = (1 - i_delta) * (1 - j_delta) * readImage[(i_index * imageWidth + j_index) * bytesPerPixel + k];
				b = (1 - i_delta) * j_delta * readImage[(i_index * imageWidth + r_index) * bytesPerPixel + k];
				c = i_delta * (1 - j_delta) * readImage[(b_index * imageWidth + j_index) * bytesPerPixel + k];
				d = i_delta * j_delta * readImage[(b_index * imageWidth + r_index) * bytesPerPixel + k];

				resizedImage[(i * resizeWidth + j) * bytesPerPixel + k] = (a + b + c + d);

				out << count++ << endl;

			}
		}
	}


	if (fopen_s(&imageWrite, argv[2], "wb"))
	{
		cout << "Cannout open output file" << endl;
		exit(1);
	}
	fwrite(&resizedImage[0], sizeof(unsigned char), resizeHeight * resizeWidth * bytesPerPixel, imageWrite);
	fclose(imageWrite);
}

void cmyk(int argc, char* argv[])
{

	int imageWidth, imageHeight, bytesPerPixel;
	int row_it, col_it, pixel_it;
	int rLoc, gLoc, bLoc, valTo;
	FILE* imageRead = 0;
	FILE* imageWrite = 0;
	FILE* imageWriteCyan = 0;
	FILE* imageWriteMagenta = 0;
	FILE* imageWriteYellow = 0;
	ofstream out("cmykOut.txt");


	if (argc < 6)
	{
		cout << argv[0] << "input_image.raw " << "output_image.raw " << "width " << "height " << "bytesPerPixel " << endl;
		exit(1);
	}

	imageWidth = atoi(argv[3]);
	imageHeight = atoi(argv[4]);
	bytesPerPixel = atoi(argv[5]);

	vector<unsigned char> rgbImage(imageHeight * imageWidth * bytesPerPixel);
	vector<unsigned char> cmykImage(imageHeight * imageWidth * bytesPerPixel);
	vector<unsigned char> cyanImage(imageHeight * imageWidth);
	vector<unsigned char> magentaImage(imageHeight * imageWidth);
	vector<unsigned char> yellowImage(imageHeight * imageWidth);

	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << argv[1] << endl;
		cout << "Cannot open Input image" << endl;
		exit(1);
	}
	fread(&rgbImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageRead);
	fclose(imageRead);
	int count = 0;
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			rLoc = (row_it * imageWidth + col_it) * bytesPerPixel + 0;
			gLoc = (row_it * imageWidth + col_it) * bytesPerPixel + 1;
			bLoc = (row_it * imageWidth + col_it) * bytesPerPixel + 2;

			cmykImage[rLoc] = (unsigned char)(255 - rgbImage[rLoc]);
			cmykImage[gLoc] = (unsigned char)(255 - rgbImage[gLoc]);
			cmykImage[bLoc] = (unsigned char)(255 - rgbImage[bLoc]);

			valTo = row_it * imageWidth + col_it;

			cyanImage[valTo] = (unsigned char)(255 - rgbImage[rLoc]);
			magentaImage[valTo] = (unsigned char)(255 - rgbImage[gLoc]);
			yellowImage[valTo] = (unsigned char)(255 - rgbImage[bLoc]);
			count++;

		}
	}

	out << count;
	// Grayscale images and final output printed

	if (fopen_s(&imageWrite, argv[2], "wb"))
	{
		cout << "Cannout open output file" << endl;
		exit(1);
	}
	fwrite(&cmykImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageWrite);
	fclose(imageWrite);

	if (fopen_s(&imageWriteCyan, "cyan.raw", "wb"))
	{
		cout << "Cannout open output file" << endl;
		exit(1);
	}
	fwrite(&cyanImage[0], sizeof(unsigned char), imageHeight * imageWidth, imageWriteCyan);
	fclose(imageWriteCyan);


	if (fopen_s(&imageWriteMagenta, "magenta.raw", "wb"))
	{
		cout << "Cannout open output file" << endl;
		exit(1);
	}
	fwrite(&magentaImage[0], sizeof(unsigned char), imageHeight * imageWidth, imageWriteMagenta);
	fclose(imageWriteMagenta);

	if (fopen_s(&imageWriteYellow, "yellow.raw", "wb"))
	{
		cout << "Cannout open output file" << endl;
		exit(1);
	}
	fwrite(&cyanImage[0], sizeof(unsigned char), imageHeight * imageWidth, imageWriteYellow);
	fclose(imageWriteYellow);
}


void hsl(int argc, char* argv[])
{

	int imageWidth, imageHeight, bytesPerPixel;
	int row_it, col_it, pixel_it;
	int rLoc, gLoc, bLoc;
	int valToH = 0, valToS = 0, valToL = 0;
	float rNorm = 0, gNorm = 0, bNorm = 0;
	FILE* imageRead = 0;
	FILE* imageWriteHue = 0;
	FILE* imageWriteLum = 0;
	FILE* imageWriteSat = 0;
	ofstream out("hslOut.txt");


	if (argc < 6)
	{
		cout << argv[0] << "input_image.raw " << "output_image.raw " << "width " << "height " << "bytesPerPixel " << endl;
		exit(1);
	}

	imageWidth = atoi(argv[3]);
	imageHeight = atoi(argv[4]);
	bytesPerPixel = atoi(argv[5]);

	vector <unsigned char> rgbImage(imageHeight * imageWidth * bytesPerPixel);
	vector <unsigned char> hslImage(imageHeight * imageWidth * bytesPerPixel); //Is this reqd?
	vector <float> hueTemp(imageHeight * imageWidth);
	vector <float> luminanceTemp(imageHeight * imageWidth);
	vector <float> saturationTemp(imageHeight * imageWidth);
	vector <unsigned char> hueImage(imageHeight * imageWidth);
	vector <unsigned char> luminanceImage(imageHeight * imageWidth);
	vector <unsigned char> saturationImage(imageHeight * imageWidth);


	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << argv[1] << endl;
		cout << "Cannot open Input image" << endl;
		exit(1);
	}
	fread(&rgbImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageRead);
	fclose(imageRead);

	float maxVal, minVal;
	float chroma, hueVal, luminenceVal, saturationVal;
	float ifMequalsR, ifMequalsG, ifMequalsB;

	out << "Hue " << "Saturation " << "Luminensce " << endl;
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			rLoc = (row_it * imageWidth + col_it) * bytesPerPixel + 0;
			gLoc = (row_it * imageWidth + col_it) * bytesPerPixel + 1;
			bLoc = (row_it * imageWidth + col_it) * bytesPerPixel + 2;

			rNorm = rgbImage[rLoc];
			gNorm = rgbImage[gLoc];
			bNorm = rgbImage[bLoc];

			// Calculate priors
			maxVal = max(max(rNorm, gNorm), bNorm);
			minVal = min(min(rNorm, gNorm), bNorm);
			chroma = maxVal - minVal;

			//For black pixels
			if (maxVal == 0.0 && minVal == 0.0)
			{
				hueVal = 0;
				luminenceVal = 0;
				saturationVal = 0;
			}

			// Calculate priors for HSV and HSV values 
			else
			{
				ifMequalsR = (gNorm - bNorm) / chroma;
				ifMequalsR = (ifMequalsR >= 6.0) ? ifMequalsR - 6.0 : ifMequalsR;
				ifMequalsR *= 60.0;
				ifMequalsG = 60.0 * (((bNorm - rNorm) / chroma) + 2.0);
				ifMequalsB = 60.0 * (((rNorm - gNorm) / chroma) + 4.0);


				hueVal = (maxVal == minVal) ? 0 : (maxVal == rNorm) ? ifMequalsR : (maxVal == gNorm) ? ifMequalsG : (maxVal == bNorm) ? ifMequalsB : 0.0;
				luminenceVal = ((maxVal + minVal) / 2.0);
				saturationVal = (luminenceVal > 0 && luminenceVal <= 0.5) ? (chroma / (2.0* luminenceVal)) : (luminenceVal == 0) ? 0.0 : (luminenceVal == 1) ? 1.00 : (chroma / (2.0 - (2.0 * luminenceVal))); // Take care if L = 1, Saturation 0 or 1?
			}

			//[FIXME] Can a HSL image be written like BGR?
			hueTemp[valToH++] = (hueVal);
			luminanceTemp[valToL++] = (luminenceVal);
			saturationTemp[valToS++] = (saturationVal);
			out << hueVal << " " << saturationVal << " " << luminenceVal << " " << valToH << "" << endl;
		}
	}

	int index = 0;
	float minHue = hueTemp[0];
	float maxHue = hueTemp[0];
	float minSat = saturationTemp[0];
	float maxSat = saturationTemp[0];
	float minLum = luminanceTemp[0];
	float maxLum = luminanceTemp[0];

	// Calc max and min values of the buffers
	while (index < (hueImage.size() - 1))
	{
		minHue = min(minHue, hueTemp[index + 1]);
		maxHue = max(maxHue, hueTemp[index + 1]);

		minSat = min(minSat, saturationTemp[index + 1]);
		maxSat = max(maxSat, saturationTemp[index + 1]);

		minLum = min(minLum, luminanceTemp[index + 1]);
		maxLum = max(maxLum, luminanceTemp[index + 1]);

		index++;
	}

	// Normalize to 0 -255 and typecast for display
	for (int index = 0; index < hueImage.size(); index++)
	{
		hueImage[index] = (unsigned char)round((hueTemp[index] - minHue) * (255.0 / (maxHue - minHue)));
		saturationImage[index] = (unsigned char)round((saturationTemp[index] - minSat) * (255.0 / (maxSat - minSat)));
		luminanceImage[index] = (unsigned char)round((luminanceTemp[index] - minLum) * (255.0 / (maxLum - minLum)));
	}

	// Write to 3 output files 
	//
	if (fopen_s(&imageWriteHue, "hue.raw", "wb"))
	{
		cout << "Cannout open output file" << endl;
		exit(1);
	}
	fwrite(&hueImage[0], sizeof(unsigned char), imageHeight * imageWidth, imageWriteHue);
	fclose(imageWriteHue);



	if (fopen_s(&imageWriteSat, "saturation.raw", "wb"))
	{
		cout << "Cannout open output file" << endl;
		exit(1);
	}
	fwrite(&saturationImage[0], sizeof(unsigned char), imageHeight * imageWidth, imageWriteSat);
	fclose(imageWriteSat);



	if (fopen_s(&imageWriteLum, "luminance.raw", "wb"))
	{
		cout << "Cannout open output file" << endl;
		exit(1);
	}
	fwrite(&luminanceImage[0], sizeof(unsigned char), imageHeight * imageWidth, imageWriteLum);
	fclose(imageWriteLum);


}
void histEqualisationGrayBucket(int argc, char* argv[])
{
	FILE* imageRead = 0;
	FILE* imageWrite = 0;
	int bytesPerPixel, imageWidth, imageHeight, imageSize;
	int row_it = 0, col_it = 0, pixel_it;
	int valTo = 0, valFrom = 0;
	int index = 0;
	double minImage, maxImage, minHist, maxHist;
	ofstream out("outHistGray.txt");

	if (argc < 6)
	{
		cerr << "Usage:" << argv[0] << "inputImage.raw" << "outputImage.raw" << "breadth" << "length" << "BytesPerPixel" << endl;
		exit(0);
	}

	imageWidth = atoi(argv[3]);
	imageHeight = atoi(argv[4]);
	bytesPerPixel = atoi(argv[5]);

	vector <unsigned char> readImage(bytesPerPixel* imageHeight* imageWidth);
	vector <unsigned char> histImage(bytesPerPixel* imageHeight* imageWidth);
	vector <double> countVal(level);  //On;ly 256 levels
	vector <double> countoVal(level);
	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << "Cannot open file" << argv[1] << endl;
		exit(1);
	}

	fread((unsigned char*)&readImage[0], sizeof(unsigned char), imageHeight* imageWidth*bytesPerPixel, imageRead);
	fclose(imageRead);



	imageSize = imageHeight * imageWidth * bytesPerPixel;

	vector <double> bucketMap(bytesPerPixel * imageHeight * imageWidth);
	vector <double> bucketNumMap(bytesPerPixel * level);
	imageSize = imageHeight * imageWidth * bytesPerPixel;


	double streched;
	double intensity;
	int bucketSize = imageHeight * imageWidth * bytesPerPixel / 256;
	int bucketNum = 0, inBucket = 0;
	double valToBucket = readImage[0];
	unsigned char temp;

	for (int i = 0; i < imageHeight * imageWidth * bytesPerPixel; i++)
	{
		bucketMap[i] = floor(i / bucketSize);
	}


	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			intensity = readImage[row_it * imageWidth + col_it];
			countVal[intensity]++;

		}
	}


	ofstream outGrayHist("plotGrayHist.txt");


	// Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outGrayHist << index << " " << countVal[index] << endl;

	}


	for (int i = 1; i < 256; i++)
	{
		countVal[i] += countVal[i - 1];
	}

	int indexColor, indexMap, count = 0;

	for (int k = 0; k < bytesPerPixel; k++)
	{
		bucketNumMap[k * bytesPerPixel + 0] = 0;
		for (int i = 1; i < level; i++)
			bucketNumMap[k * bytesPerPixel + i] = countVal[k * level + (i - 1)];

		for (row_it = 0; row_it < imageHeight; row_it++)
		{
			for (col_it = 0; col_it < imageWidth; col_it++)
			{
				count++;
				valTo = (row_it * imageWidth + col_it) * bytesPerPixel + k;
				indexColor = readImage[valTo];
				indexMap = bucketNumMap[k * bytesPerPixel + indexColor]++;
				bucketMap[indexMap] = (bucketMap[indexMap] > 255) ? 255 : round(bucketMap[indexMap]);
				histImage[valTo] = (unsigned char)bucketMap[indexMap];
				out << count << endl;
			}
		}
	}


	// Plot CDF
	ofstream outGrayCDF("plotGrayCDF.txt");

	for (int index = 0; index < 256; index++)
	{
		outGrayCDF << index << " " << countVal[index] << endl;
	}



	// Plot histogram of resulting output
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			intensity = histImage[(row_it * imageWidth + col_it) * bytesPerPixel + 0];
			countoVal[intensity]++;
		}
	}

	ofstream outoGrayHist("plotModifiedGrayHist.txt");

	// Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outoGrayHist << index << " " << countoVal[index] << endl;
	}


	if ((fopen_s(&imageWrite, argv[2], "wb")))

	{
		cout << "Cannot open output file " << endl;
	}

	fwrite((unsigned char*)&histImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageWrite);
	fclose(imageWrite);
}


void histEqualisationGrayTF(int argc, char* argv[])
{
	FILE* imageRead = 0;
	FILE* imageWrite = 0;
	int bytesPerPixel, imageWidth, imageHeight, imageSize;
	int row_it = 0, col_it = 0, pixel_it;
	int valTo = 0, valFrom = 0;
	int index = 0;
	double minImage, maxImage, minHist, maxHist;
	ofstream out("outHistGray.txt");

	if (argc < 6)
	{
		cerr << "Usage:" << argv[0] << "inputImage.raw" << "outputImage.raw" << "breadth" << "length" << "BytesPerPixel" << endl;
		exit(0);
	}

	imageWidth = atoi(argv[3]);
	imageHeight = atoi(argv[4]);
	bytesPerPixel = atoi(argv[5]);

	vector <unsigned char> readImage(bytesPerPixel* imageHeight* imageWidth);
	vector <unsigned char> histImage(bytesPerPixel* imageHeight* imageWidth);
	vector <double> countVal(level);  //On;ly 256 levels
	vector <double> countoVal(level);

	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << "Cannot open file" << argv[1] << endl;
		exit(1);
	}

	fread((unsigned char*)&readImage[0], sizeof(unsigned char), imageHeight* imageWidth*bytesPerPixel, imageRead);
	fclose(imageRead);


	maxHist = 255;
	minHist = 0;
	double streched;
	double intensity;

	// Count intensities
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			intensity = readImage[row_it * imageWidth + col_it];
			countVal[intensity]++;
		}
	}



	ofstream outGrayHist("plotGrayHist.txt");


	// Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outGrayHist << index << " " << countVal[index] << endl;

	}


	// Normalize histogram by the image size

	for (row_it = 0; row_it < 256; row_it++)
	{
		countVal[row_it] /= imageHeight * imageWidth;
	}

	//Caclculate the cdf
	for (int index = 1; index < 256; index++)
	{
		countVal[index] = (countVal[index] == 1) ? 1 : (countVal[index - 1] + countVal[index]);
	}

	// Normalize to the range
	for (int index = 0; index < 256; index++)
	{
		countVal[index] *= 255;

		// Clip the values
		countVal[index] = countVal[index] > 255 ? 255 : countVal[index];
	}


	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			valTo = row_it * imageWidth + col_it;
			valFrom = valTo;
			histImage[valTo] = (unsigned char)countVal[readImage[valFrom]];

		}
	}



	// Plot CDF
	ofstream outGrayCDF("plotGrayCDF.txt");

	for (int index = 0; index < 256; index++)
	{
		outGrayCDF << index << " " << countVal[index] << endl;
	}



	// Plot histogram of resulting output
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			intensity = histImage[(row_it * imageWidth + col_it) * bytesPerPixel + 0];
			countoVal[intensity]++;
		}
	}

	ofstream outoGrayHist("plotModifiedGrayHist.txt");

	// Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outoGrayHist << index << " " << countoVal[index] << endl;
	}






	if ((fopen_s(&imageWrite, argv[2], "wb")))
	{
		cout << "Cannot open output file " << endl;
	}

	fwrite((unsigned char*)&histImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageWrite);
	fclose(imageWrite);
}

void histEqualisationColorBucket(int argc, char* argv[])
{
	FILE* imageRead = 0;
	FILE* imageWrite = 0;
	int bytesPerPixel, imageWidth, imageHeight, imageSize;
	int row_it = 0, col_it = 0, pixel_it;
	int valTo = 0, valFrom = 0;
	int index = 0;
	double minImage, maxImage, minHist, maxHist;
	ofstream out("outHistGray.txt");

	if (argc < 6)
	{
		cerr << "Usage:" << argv[0] << "inputImage.raw" << "outputImage.raw" << "breadth" << "length" << "BytesPerPixel" << endl;
		exit(0);
	}

	imageWidth = atoi(argv[3]);
	imageHeight = atoi(argv[4]);
	bytesPerPixel = atoi(argv[5]);

	vector <unsigned char> readImage(bytesPerPixel* imageHeight* imageWidth);
	vector <unsigned char> histImage(bytesPerPixel* imageHeight* imageWidth);
	vector <double> countRVal(level);  //Only 256 levels
	vector <double> countGVal(level);
	vector <double> countBVal(level);
	vector <double> countoRVal(level);
	vector <double> countoGVal(level);
	vector <double> countoBVal(level);

	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << "Cannot open file" << argv[1] << endl;
		exit(1);
	}

	fread((unsigned char*)&readImage[0], sizeof(unsigned char), imageHeight* imageWidth*bytesPerPixel, imageRead);
	fclose(imageRead);

	imageSize = imageHeight * imageWidth * bytesPerPixel;

	vector <double> bucketMap(bytesPerPixel * imageHeight * imageWidth);
	vector <double> bucketNumMapR(level);
	vector <double> bucketNumMapG(level);
	vector <double> bucketNumMapB(level);
	imageSize = imageHeight * imageWidth * bytesPerPixel;


	double streched;
	double rIntensity, gIntensity, bIntensity;
	int bucketSize = imageHeight * imageWidth / 256;
	int bucketNum = 0, inBucket = 0;
	double valToBucket = readImage[0];
	unsigned char temp;

	for (int i = 0; i < imageHeight * imageWidth * bytesPerPixel; i++)
	{
		bucketMap[i] = floor(i / bucketSize);
	}


	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			valTo = row_it * imageWidth + col_it;
			rIntensity = readImage[valTo * bytesPerPixel + 0];
			gIntensity = readImage[valTo * bytesPerPixel + 1];
			bIntensity = readImage[valTo * bytesPerPixel + 2];
			countRVal[rIntensity]++;
			countGVal[gIntensity]++;
			countBVal[bIntensity]++;

		}
	}

	ofstream outRHist("plotRHist.txt");
	ofstream outGHist("plotGHist.txt");
	ofstream outBHist("plotBHist.txt");

	// Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outRHist << index << " " << countRVal[index] << endl;
		outGHist << index << " " << countGVal[index] << endl;
		outBHist << index << " " << countBVal[index] << endl;
	}




	for (int i = 1; i < 256; i++)
	{
		countRVal[i] += countRVal[i - 1];
		countGVal[i] += countGVal[i - 1];
		countBVal[i] += countBVal[i - 1];

	}


	int indexColorR, indexColorG, indexColorB, indexMapR, indexMapG, indexMapB, count = 0;


	bucketNumMapR[0] = 0;
	bucketNumMapG[0] = 0;
	bucketNumMapB[0] = 0;

	for (int i = 1; i < level; i++)
	{
		bucketNumMapR[(i)] = countRVal[(i - 1)];
		bucketNumMapG[(i)] = countGVal[(i - 1)];
		bucketNumMapB[(i)] = countBVal[(i - 1)];

		out << "R" << " " << (i) << " " << bucketNumMapR[(i)] << endl;
		out << "G" << " " << (i) << " " << bucketNumMapG[(i)] << endl;
		out << "B" << " " << (i) << " " << bucketNumMapR[(i)] << endl;

	}
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			count++;
			valTo = (row_it * imageWidth + col_it);
			indexColorR = readImage[valTo * bytesPerPixel + 0];
			indexColorG = readImage[valTo * bytesPerPixel + 1];
			indexColorB = readImage[valTo * bytesPerPixel + 2];
			indexMapR = bucketNumMapR[indexColorR]++;
			indexMapG = bucketNumMapG[indexColorG]++;
			indexMapB = bucketNumMapB[indexColorB]++;
			(bucketMap[indexMapR]) = (bucketMap[indexMapR] > 255) ? 255 : round(bucketMap[indexMapR]);
			(bucketMap[indexMapG]) = (bucketMap[indexMapG] > 255) ? 255 : round(bucketMap[indexMapG]);
			(bucketMap[indexMapB]) = (bucketMap[indexMapB] > 255) ? 255 : round(bucketMap[indexMapB]);
			histImage[valTo * bytesPerPixel + 0] = (unsigned char)bucketMap[indexMapR];
			histImage[valTo * bytesPerPixel + 1] = (unsigned char)bucketMap[indexMapG];
			histImage[valTo * bytesPerPixel + 2] = (unsigned char)bucketMap[indexMapB];
			//out << count << endl;
		}
	}


	//Plot CDF
	ofstream outRCDF("plotRCDF.txt");
	ofstream outGCDF("plotGCDF.txt");
	ofstream outBCDF("plotBCDF.txt");

	for (int index = 0; index < 256; index++)
	{

		// Write to output
		outRCDF << index << " " << bucketNumMapR[index] << endl;
		outGCDF << index << " " << bucketNumMapG[index] << endl;
		outBCDF << index << " " << bucketNumMapB[index] << endl;
	}



	// Plot histogram of resulting output
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{

			rIntensity = histImage[(row_it * imageWidth + col_it) * bytesPerPixel + 0];
			countoRVal[rIntensity]++;

			gIntensity = histImage[(row_it * imageWidth + col_it) * bytesPerPixel + 1];
			countoGVal[gIntensity]++;

			bIntensity = histImage[(row_it * imageWidth + col_it) * bytesPerPixel + 2];
			countoBVal[bIntensity]++;

		}
	}

	ofstream outoRHist("plotModifiedRHist.txt");
	ofstream outoGHist("plotModifiedGHist.txt");
	ofstream outoBHist("plotModifiedBHist.txt");

	// Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outoRHist << index << " " << countoRVal[index] << endl;
		outoGHist << index << " " << countoGVal[index] << endl;
		outoBHist << index << " " << countoBVal[index] << endl;
	}



	if ((fopen_s(&imageWrite, argv[2], "wb")))

	{
		cout << "Cannot open output file " << endl;
	}

	fwrite((unsigned char*)&histImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageWrite);
	fclose(imageWrite);
}

void histEqualisationColorTF(int argc, char* argv[])
{
	FILE* imageRead = 0;
	FILE* imageWrite = 0;
	int bytesPerPixel, imageWidth, imageHeight, imageSize;
	int row_it = 0, col_it = 0, pixel_it;
	int valToR = 0, valToG = 0, valToB = 0, valFromR = 0, valFromG = 0, valFromB = 0;
	int index = 0;
	double minImage, maxImage, minHist, maxHist;
	ofstream out("outHistGray.txt");

	if (argc < 6)
	{
		cerr << "Usage:" << argv[0] << "inputImage.raw" << "outputImage.raw" << "breadth" << "length" << "BytesPerPixel" << endl;
		exit(0);
	}

	imageWidth = atoi(argv[3]);
	imageHeight = atoi(argv[4]);
	bytesPerPixel = atoi(argv[5]);

	vector <unsigned char> readImage(bytesPerPixel* imageHeight* imageWidth);
	vector <unsigned char> histImage(bytesPerPixel* imageHeight* imageWidth);
	vector <double> countRVal(level);
	vector <double> countGVal(level);
	vector <double> countBVal(level);
	vector <double> countoRVal(level);
	vector <double> countoGVal(level);
	vector <double> countoBVal(level);


	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << "Cannot open file" << argv[1] << endl;
		exit(1);
	}

	fread((unsigned char*)&readImage[0], sizeof(unsigned char), imageHeight* imageWidth*bytesPerPixel, imageRead);
	fclose(imageRead);


	double rIntensity, gIntensity, bIntensity;

	// Count intensities
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{

			rIntensity = readImage[(row_it * imageWidth + col_it) * bytesPerPixel + 0];
			countRVal[rIntensity]++;

			gIntensity = readImage[(row_it * imageWidth + col_it) * bytesPerPixel + 1];
			countGVal[gIntensity]++;

			bIntensity = readImage[(row_it * imageWidth + col_it) * bytesPerPixel + 2];
			countBVal[bIntensity]++;

		}
	}

	ofstream outRHist("plotRHist.txt");
	ofstream outGHist("plotGHist.txt");
	ofstream outBHist("plotBHist.txt");

	// Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outRHist << index << " " << countRVal[index] << endl;
		outGHist << index << " " << countGVal[index] << endl;
		outBHist << index << " " << countBVal[index] << endl;
	}


	// Normalize histogram by the image size

	for (row_it = 0; row_it < 256; row_it++)
	{
		countRVal[row_it] /= imageHeight * imageWidth;
		countGVal[row_it] /= imageHeight * imageWidth;
		countBVal[row_it] /= imageHeight * imageWidth;
	}

	//Caclculate the cdf
	for (int index = 1; index < 256; index++)
	{
		countRVal[index] = (countRVal[index] == 1) ? 1 : (countRVal[index - 1] + countRVal[index]);
		countGVal[index] = (countGVal[index] == 1) ? 1 : (countGVal[index - 1] + countGVal[index]);
		countBVal[index] = (countBVal[index] == 1) ? 1 : (countBVal[index - 1] + countBVal[index]);
	}




	// Normalize to the range
	for (int index = 0; index < 256; index++)
	{
		countRVal[index] *= 255;
		countGVal[index] *= 255;
		countBVal[index] *= 255;
		// Clip the values
		countRVal[index] = countRVal[index] > 255 ? 255 : countRVal[index];
		countGVal[index] = countGVal[index] > 255 ? 255 : countGVal[index];
		countBVal[index] = countBVal[index] > 255 ? 255 : countBVal[index];
	}

	//Plot CDF



	ofstream outRCDF("plotRCDF.txt");
	ofstream outGCDF("plotGCDF.txt");
	ofstream outBCDF("plotBCDF.txt");

	for (int index = 0; index < 256; index++)
	{

		//[MAYDAY : Apparently not need this ?: Clarify with TA]

		/*// Normalize the CDF over the image
		countRVal[index] /= imageHeight * imageWidth;
		countGVal[index] /= imageHeight * imageWidth;
		countBVal[index] /= imageHeight * imageWidth;

		// Clip the CDF between 0 and 1

		countRVal[index] = (countRVal[index] > 255) ? 255 : (countRVal[index] < 0) ? 0 : countRVal[index];
		countGVal[index] = (countGVal[index] > 255) ? 255 : (countGVal[index] < 0) ? 0 : countGVal[index];
		countBVal[index] = (countRVal[index] > 255) ? 255 : (countBVal[index] < 0) ? 0 : countBVal[index];*/

		// Write to output
		outRCDF << index << " " << countRVal[index] << endl;
		outGCDF << index << " " << countGVal[index] << endl;
		outBCDF << index << " " << countBVal[index] << endl;
	}


	// Write to data to output buffer
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			valToR = (row_it * imageWidth + col_it) * bytesPerPixel + 0;
			valToG = (row_it * imageWidth + col_it) * bytesPerPixel + 1;
			valToB = (row_it * imageWidth + col_it) * bytesPerPixel + 2;
			valFromR = valToR;
			valFromG = valToG;
			valFromB = valToB;
			histImage[valToR] = (unsigned char)countRVal[readImage[valFromR]];
			histImage[valToG] = (unsigned char)countGVal[readImage[valFromG]];
			histImage[valToB] = (unsigned char)countBVal[readImage[valFromB]];
		}
	}

	// Plot histogram of resulting output
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{

			rIntensity = histImage[(row_it * imageWidth + col_it) * bytesPerPixel + 0];
			countoRVal[rIntensity]++;

			gIntensity = histImage[(row_it * imageWidth + col_it) * bytesPerPixel + 1];
			countoGVal[gIntensity]++;

			bIntensity = histImage[(row_it * imageWidth + col_it) * bytesPerPixel + 2];
			countoBVal[bIntensity]++;

		}
	}

	ofstream outoRHist("plotModifiedRHist.txt");
	ofstream outoGHist("plotModifiedGHist.txt");
	ofstream outoBHist("plotModifiedBHist.txt");

	// Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outoRHist << index << " " << countoRVal[index] << endl;
		outoGHist << index << " " << countoGVal[index] << endl;
		outoBHist << index << " " << countoBVal[index] << endl;
	}

	if ((fopen_s(&imageWrite, argv[2], "wb")))
	{
		cout << "Cannot open output file " << endl;
	}

	fwrite((unsigned char*)&histImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageWrite);
	fclose(imageWrite);

}

void specialEffect(int argc, char* argv[])
{
	int imageWidth, imageHeight, bytesPerPixel;
	int row_it, col_it, pixel_it;
	int valTo, valFrom;
	FILE* imageRead = 0;
	FILE* imageWrite = 0;
	ofstream out("filmOut.txt");
	vector<matchPoints>sourceR;
	vector<matchPoints>sourceG;
	vector<matchPoints>sourceB;
	vector<matchPoints>desiredR;
	vector<matchPoints>desiredG;
	vector<matchPoints>desiredB;
	vector<matchPoints>afterMatchR;
	vector<matchPoints>afterMatchG;
	vector<matchPoints>afterMatchB;
	vector <double> matchCDF(level);
	vector <double> countGaussian(level);
	vector<double> desiredHist(level);

	if (argc < 6)
	{
		cout << argv[0] << "input_image.raw " << "output_image.raw " << "width " << "height " << "bytesPerPixel " << endl;
		exit(1);
	}

	imageWidth = atoi(argv[3]);
	imageHeight = atoi(argv[4]);
	bytesPerPixel = atoi(argv[5]);


	vector<unsigned char> readImage(imageHeight * imageWidth * bytesPerPixel);
	vector<unsigned char> effectImage(imageHeight * imageWidth * bytesPerPixel);
	vector <double> countRVal(level);
	vector <double> countGVal(level);
	vector <double> countBVal(level);

	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << argv[1] << endl;
		cout << "Cannot open Input image" << endl;
		exit(1);
	}
	fread(&readImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageRead);
	fclose(imageRead);

	unsigned char intensityR, intensityG, intensityB;
	double totalIntensity;
	int imageSize = imageHeight * imageWidth * bytesPerPixel;


	for (int index = 0; index < 256; index++)
	{
		sourceR.push_back(matchPoints());
		desiredR.push_back(matchPoints());
		afterMatchR.push_back(matchPoints());


		sourceG.push_back(matchPoints());
		desiredG.push_back(matchPoints());
		afterMatchG.push_back(matchPoints());


		sourceB.push_back(matchPoints());
		desiredB.push_back(matchPoints());
		afterMatchB.push_back(matchPoints());
	}


	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			intensityR = readImage[(row_it * imageWidth + col_it) + bytesPerPixel + 0];
			intensityG = readImage[(row_it * imageWidth + col_it) + bytesPerPixel + 1];
			intensityB = readImage[(row_it * imageWidth + col_it) + bytesPerPixel + 2];
			countRVal[intensityR]++;
			countRVal[intensityG]++;
			countRVal[intensityB]++;
		}
	}


	ofstream outRHist("plotRHist.txt");
	ofstream outGHist("plotGHist.txt");
	ofstream outBHist("plotBHist.txt");


	// Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outRHist << index << " " << countRVal[index] << endl;
		outGHist << index << " " << countGVal[index] << endl;
		outBHist << index << " " << countBVal[index] << endl;

	}

	// Normalize histogram by the image size
	for (int index = 0; index < 256; index++)
	{
		countRVal[index] /= imageWidth * imageHeight;
		countGVal[index] /= imageWidth * imageHeight;
		countBVal[index] /= imageWidth * imageHeight;
	}


	//Caclculate the cdf
	for (int index = 1; index < 256; index++)
	{
		countRVal[index] = (countRVal[index] == 1) ? 1 : (countRVal[index - 1] + countRVal[index]);
		countGVal[index] = (countGVal[index] == 1) ? 1 : (countGVal[index - 1] + countGVal[index]);
		countBVal[index] = (countBVal[index] == 1) ? 1 : (countBVal[index - 1] + countBVal[index]);
	}

	// Normalize to the range
	for (int index = 0; index < 256; index++)
	{
		countRVal[index] *= 255;
		countGVal[index] *= 255;
		countBVal[index] *= 255;

		// Clip the values
		countRVal[index] = countRVal[index] > 255 ? 255 : countRVal[index];
		countGVal[index] = countGVal[index] > 255 ? 255 : countGVal[index];
		countBVal[index] = countBVal[index] > 255 ? 255 : countBVal[index];

		sourceR[index].pixelIntensity = index;
		sourceR[index].cdf = countRVal[index];

		sourceG[index].pixelIntensity = index;
		sourceG[index].cdf = countGVal[index];

		sourceB[index].pixelIntensity = index;
		sourceB[index].cdf = countBVal[index];
	}

	ifstream inR("plotDRCDF.txt");
	ifstream inG("plotDGCDF.txt");
	ifstream inB("plotDBCDF.txt");

	for (int index = 0; index < level; index++)
	{

		inR >> desiredR[index].pixelIntensity >> desiredR[index].cdf;
		inG >> desiredG[index].pixelIntensity >> desiredG[index].cdf;
		inB >> desiredB[index].pixelIntensity >> desiredB[index].cdf;

	}

	vector<int> countAfterHistMatch(level);

	for (int index = 0; index < level; index++)
	{
		if ((desiredR[index].cdf - sourceR[index].cdf) < 0)
		{
			afterMatchR[index].cdf = int(sourceR[index].cdf);
			afterMatchR[index].pixelIntensity = desiredR[int(afterMatchR[index].cdf)].pixelIntensity;
		}
		else
		{
			afterMatchR[index].cdf = desiredR[index].cdf;
			afterMatchR[index].pixelIntensity = desiredR[index].pixelIntensity;
		}


		if ((desiredG[index].cdf - sourceG[index].cdf) < 0)
		{
			afterMatchG[index].cdf = int(sourceG[index].cdf);
			afterMatchG[index].pixelIntensity = desiredG[int(afterMatchG[index].cdf)].pixelIntensity;
		}
		else
		{
			afterMatchG[index].cdf = desiredG[index].cdf;
			afterMatchG[index].pixelIntensity = desiredG[index].pixelIntensity;
		}


		if ((desiredB[index].cdf - sourceB[index].cdf) < 0)
		{
			afterMatchB[index].cdf = int(sourceB[index].cdf);
			afterMatchB[index].pixelIntensity = desiredB[int(afterMatchB[index].cdf)].pixelIntensity;
		}
		else
		{
			afterMatchB[index].cdf = desiredB[index].cdf;
			afterMatchB[index].pixelIntensity = desiredB[index].pixelIntensity;
		}

	}



	for (int row_it = 0; row_it < imageHeight; row_it++)
	{
		for (int col_it = 0; col_it < imageWidth; col_it++)
		{
			effectImage[(row_it * imageWidth + col_it) * bytesPerPixel + 0] = afterMatchR[readImage[(row_it * imageWidth + col_it) * bytesPerPixel + 0]].pixelIntensity;
			effectImage[(row_it * imageWidth + col_it) * bytesPerPixel + 1] = afterMatchG[readImage[(row_it * imageWidth + col_it) * bytesPerPixel + 1]].pixelIntensity;
			effectImage[(row_it * imageWidth + col_it) * bytesPerPixel + 2] = afterMatchB[readImage[(row_it * imageWidth + col_it) * bytesPerPixel + 2]].pixelIntensity;
		}
	}



	vector<double>countoRVal(imageHeight * imageWidth);
	vector<double>countoGVal(imageHeight * imageWidth);
	vector<double>countoBVal(imageHeight * imageWidth);

	for (int row_it = 0; row_it < imageHeight; row_it++)
	{
		for (int col_it = 0; col_it < imageWidth; col_it++)
		{
			intensityR = effectImage[(row_it * imageWidth + col_it) * bytesPerPixel + 0];
			intensityG = effectImage[(row_it * imageWidth + col_it) * bytesPerPixel + 1];
			intensityB = effectImage[(row_it * imageWidth + col_it) * bytesPerPixel + 2];
			countoRVal[intensityR]++;
			countoGVal[intensityG]++;
			countoBVal[intensityB]++;
		}
	}

	ofstream outoRHist("plotModifiedRHist.txt");
	ofstream outoGHist("plotModifiedGHist.txt");
	ofstream outoBHist("plotModifiedBHist.txt");

	//[FIXME] Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outoRHist << index << " " << countoRVal[index] << endl;
		outoGHist << index << " " << countoGVal[index] << endl;
		outoBHist << index << " " << countoBVal[index] << endl;
	}



	// Open out file for writing 
	if (fopen_s(&imageWrite, argv[2], "wb"))
	{
		cout << "Cannout open output file" << endl;
		exit(1);
	}
	fwrite(&effectImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageWrite);
	fclose(imageWrite);

}


void histTransform(int argc, char* argv[])
{
	FILE* imageRead = 0;
	FILE* imageWrite = 0;
	int bytesPerPixel, imageWidth, imageHeight, imageSize;
	int row_it = 0, col_it = 0, pixel_it;
	int valTo = 0, valFrom = 0;
	int index = 0;
	double minImage, maxImage, minHist, maxHist;
	ofstream out("outHistGray.txt");

	if (argc < 6)
	{
		cerr << "Usage:" << argv[0] << "inputImage.raw" << "outputImage.raw" << "breadth" << "length" << "BytesPerPixel" << endl;
		exit(0);
	}

	imageWidth = atoi(argv[3]);
	imageHeight = atoi(argv[4]);
	bytesPerPixel = atoi(argv[5]);

	vector <unsigned char> readImage(bytesPerPixel* imageHeight* imageWidth);
	vector <unsigned char> histImage(bytesPerPixel* imageHeight* imageWidth);
	vector <double> countVal(level);  //On;ly 256 levels
	vector <double> countoVal(level);
	vector<matchPoints>source;
	vector<matchPoints>gauss;
	vector<matchPoints>afterMatch;
	vector <double> matchCDF(level);
	vector <double> countGaussian(level);
	vector<double> desiredHist(level);

	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << "Cannot open file" << argv[1] << endl;
		exit(1);
	}

	fread((unsigned char*)&readImage[0], sizeof(unsigned char), imageHeight* imageWidth*bytesPerPixel, imageRead);
	fclose(imageRead);

	for (int index = 0; index < 256; index++)
	{
		source.push_back(matchPoints());
		gauss.push_back(matchPoints());
		afterMatch.push_back(matchPoints());

	}
	maxHist = 255;
	minHist = 0;
	double streched;
	double intensity;

	// Count intensities
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			intensity = readImage[row_it * imageWidth + col_it];
			countVal[intensity]++;
		}
	}


	ofstream outGrayHist("plotGrayHist.txt");


	// Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outGrayHist << index << " " << countVal[index] << endl;

	}

	// Normalize histogram by the image size
	for (int index = 0; index < 256; index++)
	{
		countVal[index] /= imageWidth * imageHeight;
	}


	//Caclculate the cdf
	for (int index = 1; index < 256; index++)
	{
		countVal[index] = (countVal[index] == 1) ? 1 : (countVal[index - 1] + countVal[index]);
	}

	// Normalize to the range
	for (int index = 0; index < 256; index++)
	{
		countVal[index] *= 255;

		// Clip the values
		countVal[index] = countVal[index] > 255 ? 255 : countVal[index];
		source[index].pixelIntensity = index;
		source[index].cdf = countVal[index];
	}



	const int noIter = 10000;
	const int stars = 100;

	default_random_engine generateNormal;
	normal_distribution<double> distribution(150.0, 30.0);

	for (int index = 0; index < noIter; index++)
	{
		double gaussVal = distribution(generateNormal);
		if ((gaussVal >= 0) && (gaussVal < level))
		{
			countGaussian[int(gaussVal)]++;

		}
	}

	for (int index = 0; index < level; index++)
	{
		gauss[index].pixelIntensity = index;
		gauss[index].cdf = countGaussian[index];
	}

	vector<int> countAfterHistMatch(level);

	for (int index = 0; index < level; index++)
	{
		if ((gauss[index].cdf - source[index].cdf) < 0)
		{
			afterMatch[index].cdf = int(source[index].cdf);
			afterMatch[index].pixelIntensity = gauss[int(afterMatch[index].cdf)].pixelIntensity;
		}
		else
		{
			afterMatch[index].cdf = gauss[index].cdf;
			afterMatch[index].pixelIntensity = gauss[index].pixelIntensity;
		}
	}



	for (int row_it = 0; row_it < imageHeight; row_it++)
	{
		for (int col_it = 0; col_it < imageWidth; col_it++)
		{
			histImage[row_it * imageWidth + col_it] = afterMatch[readImage[row_it * imageWidth + col_it]].pixelIntensity;
		}
	}


	/*for (row_it = 0; row_it < imageHeight; row_it++)
	{
	for (col_it = 0; col_it < imageWidth; col_it++)
	{
	valTo = row_it * imageWidth + col_it;
	valFrom = valTo;
	histImage[valTo] = (unsigned char)countVal[readImage[valFrom]];

	}
	}*/

	// Plot CDF
	ofstream outGrayCDF("plotGrayCDF.txt");

	for (int index = 0; index < 256; index++)
	{
		outGrayCDF << index << " " << countGaussian[index] << endl;
	}

	// Plot histogram of resulting output
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			intensity = histImage[(row_it * imageWidth + col_it) * bytesPerPixel + 0];
			countoVal[intensity]++;
		}
	}

	ofstream outoGrayHist("plotModifiedGrayHist.txt");

	// Plot the histogram of the three channels before taking the CDF

	for (int index = 0; index < 256; index++)
	{
		outoGrayHist << index << " " << countoVal[index] << endl;
	}


	if ((fopen_s(&imageWrite, argv[2], "wb")))
	{
		cout << "Cannot open output file " << endl;
	}

	fwrite((unsigned char*)&histImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageWrite);
	fclose(imageWrite);
}


unsigned char medianOf(unsigned char find[], int len)
{
	unsigned char temp;
	for (int i = 0; i < len - 1; i++)
		for (int j = i + 1; j < len; j++)
			if (find[i] > find[j]) {
				temp = find[i];
				find[i] = find[j];
				find[j] = temp;
			}
	return find[len / 2];
}



void mixNoise(int argc, char* argv[])

{
	FILE* imageRead = 0;
	FILE* imageReadOriginal = 0;
	FILE* imageWrite = 0;
	int bytesPerPixel, imageWidth, imageHeight, imageSize;
	int row_it = 0, col_it = 0, pixel_it;
	int valToR = 0, valToG = 0, valToB = 0, valFromR = 0, valFromG = 0, valFromB = 0;
	int index = 0;
	int filterSize = 0;
	double minImage, maxImage, minHist, maxHist;
	ofstream out("outHistGray.txt");

	if (argc < 8)
	{
		cerr << "Usage:" << argv[0] << "inputImage.raw" << "InputNoNoiseImg" << "outputImage.raw " << "breadth " << "length " << "BytesPerPixel " << "Filter size " << endl;
		exit(0);
	}

	imageWidth = atoi(argv[4]);
	imageHeight = atoi(argv[5]);
	bytesPerPixel = atoi(argv[6]);

	if (argc >= 6) filterSize = atoi(argv[7]);
	else filterSize = 3;

	vector <unsigned char> readImage(bytesPerPixel* imageHeight* imageWidth);
	vector <unsigned char> readImageOriginal(bytesPerPixel* imageHeight* imageWidth);
	vector <unsigned char> paddedImage(bytesPerPixel* imageHeight* imageWidth);
	vector <unsigned char> correctedImage(bytesPerPixel* imageHeight* imageWidth);



	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << "Cannot open file" << argv[1] << endl;
		exit(1);
	}

	fread((unsigned char*)&readImage[0], sizeof(unsigned char), imageHeight* imageWidth*bytesPerPixel, imageRead);
	fclose(imageRead);


	if ((fopen_s(&imageReadOriginal, argv[2], "rb")))
	{
		cout << "Cannot open file" << argv[2] << endl;
		exit(1);
	}

	fread((unsigned char*)&readImageOriginal[0], sizeof(unsigned char), imageHeight* imageWidth*bytesPerPixel, imageRead);
	fclose(imageRead);

	//gaussianFilter(&readImage[0], imageWidth, imageHeight, bytesPerPixel, filterSize);


	ofstream out1("noise.txt");
	vector<unsigned char> tempImage(imageWidth * imageHeight * bytesPerPixel);

	// Set up the gaussian window
	int centerPixel = (filterSize) / 2;
	vector<unsigned char>windowC(filterSize * filterSize);
	vector<double>gaussianWindow(filterSize * filterSize);
	double sumWindow = 0.0;
	double calcXY, calcExp;
	double sigma = 1;

	// Calculate the window values


	double convolvePixels = 0.0;
	double windowWeight = 0.0;
	double calc = 0;
	int count = 0;
	// Convole the filter to image



	// Set up the gaussian window
	int centerPixel = (filterSize) / 2;


	// Calculate the window values

	for (int row_it = 0; row_it < filterSize; row_it++)
	{
		for (int col_it = 0; col_it < filterSize; col_it++)
		{
			calcXY = (row_it - centerPixel) * (row_it - centerPixel) + (col_it - centerPixel) * (col_it - centerPixel);
			calcExp = exp(-calcXY / (2 * sigma * sigma));
			gaussianWindow[row_it * filterSize + col_it] = calcExp / (2 * M_PI * sigma);
			sumWindow += gaussianWindow[row_it * filterSize + col_it];
		}
	}

	// Normalize 

	for (int row_it = 0; row_it < filterSize; row_it++)
	{
		for (int col_it = 0; col_it < filterSize; col_it++)
		{
			gaussianWindow[row_it * filterSize + col_it] /= sumWindow;
		}
	}


	double convolvePixels = 0.0;
	double windowWeight = 0.0;
	// Convole the filter to image

	for (int row_it = 0; row_it < imageHeight; row_it++)
	{
		for (int col_it = 0; col_it < imageWidth; col_it++)
		{
			// Check for boundary conditions
			int corners[4];
			corners[0] = (row_it - centerPixel) < 0 ? 0 : (row_it - centerPixel);
			corners[1] = (col_it - centerPixel) < 0 ? 0 : (col_it - centerPixel);
			corners[2] = (row_it + centerPixel) > (imageHeight - 1) ? imageHeight : (row_it + centerPixel);
			corners[3] = (col_it + centerPixel) > (imageWidth - 1) ? imageWidth : (col_it + centerPixel);
			out << corners[0] << " " << corners[2] << " " << corners[1] << " " << corners[3] << endl;
			for (int pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
			{
				sumWindow = 0.0;
				convolvePixels = 0.0;
				windowWeight = 0.0;
				for (int filrow_it = corners[0]; filrow_it < corners[2]; filrow_it++)
				{
					for (int filcol_it = corners[1]; filcol_it < corners[3]; filcol_it++)
					{
						count++;
						out << (centerPixel + filrow_it - row_it) << " " << filterSize << " " << (centerPixel + filcol_it - col_it) << " " << "Gaussian start" << endl;
						windowWeight = gaussianWindow[(centerPixel + filrow_it - row_it) * filterSize + (centerPixel + filcol_it - col_it)];
						out << "Gaussian End" << endl;
						sumWindow += windowWeight;
						out << (filrow_it * imageWidth + filcol_it) << " " << bytesPerPixel << "  " << pixel_it << " " << "ReadImg Start" << endl;
						convolvePixels += windowWeight * readImage[(filrow_it * imageWidth + filcol_it) * bytesPerPixel + pixel_it];
						out << "ReadImgEnd" << endl;
						out << count << endl;
					}
				}

				tempImage[(row_it * imageWidth + col_it) * bytesPerPixel + pixel_it] = (unsigned char)round(convolvePixels / sumWindow);


			}
		}
	}













	// Median
	for (int row_it = 0; row_it < imageHeight; row_it++)
	{
		for (int col_it = 0; col_it < imageWidth; col_it++)
		{
			for (int pixel_it = 0; pixel_it < bytesPerPixel; pixel_it++)
			{
				count = 0;
				for (int filrow_it = row_it - filterSize / 2; filrow_it < row_it + filterSize / 2; filrow_it++)
				{
					for (int filcol_it = col_it - filterSize / 2; filcol_it < col_it + filterSize / 2; filcol_it++)
					{
						if (filrow_it >= 0 && filrow_it < imageHeight && filcol_it >= 0 && filcol_it < imageWidth)
						{
							windowC[count] = readImage[(filrow_it * (imageWidth)+filcol_it) * bytesPerPixel + pixel_it];
							count++;
						}
					}
				}
				tempImage[(row_it * imageWidth + col_it) * bytesPerPixel + pixel_it] = medianOf(&windowC[0], count);
			}
		}
	}




	double MSE[3] = { 0,0,0 };
	double PSNR[3] = { 0,0,0 };
	int newValue, oriValue;


	for (int k = 0; k < 3; k++)
	{
		for (int i = 0; i < imageHeight; i++)
		{
			for (int j = 0; j < imageWidth; j++)
			{

				newValue = tempImage[(i * imageWidth + j) * 3 + k];
				oriValue = readImageOriginal[(i * imageWidth + j) * 3 + k];

				MSE[k] += (double)(newValue - oriValue) * (newValue - oriValue) / (imageHeight * imageWidth);
			}
		}
	}
	for (int k = 0; k < 3; k++)
	{
		PSNR[k] = 10 * log10((double)(255 * 255) / MSE[k]);
		cout << "MSE = " << MSE[k] << "PSNR =" << PSNR[k] << endl;
		out << k << " " << MSE[k] << " " << PSNR[k] << endl;
	}
	printf("******\n");



	if ((fopen_s(&imageWrite, argv[3], "wb")))
	{
		cout << "Cannot open output file " << endl;
	}

	fwrite((unsigned char*)&tempImage[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageWrite);
	fclose(imageWrite);

}


/*------------------------------------3b----------------------------------------------*/
void nonLocalMean(int argc, char* argv[])
{
	FILE* imageRead = 0;
	FILE* imageWrite = 0;
	FILE* imageReadOriginal = 0;
	int bytesPerPixel, imageWidth, imageHeight, imageSize;
	int row_it = 0, col_it = 0, pixel_it;
	int valFrom = 0, valTo = 0;
	//int valToR = 0, valToG = 0, valToB = 0, valFromR = 0, valFromG = 0, valFromB = 0;
	int index = 0;
	double minImage, maxImage, minHist, maxHist;
	ofstream out("outHistGray.txt");

	if (argc < 8)
	{
		cerr << "Usage:" << argv[0] << "inputnoisyImage.raw" << "inputNoNoiseImage.raw" << "outputImage.raw" << "breadth" << "length" << "BytesPerPixel" << "Filter size" << endl;
		exit(0);
	}

	imageWidth = atoi(argv[4]);
	imageHeight = atoi(argv[5]);
	bytesPerPixel = atoi(argv[6]);
	int N;
	N = atoi(argv[7]);

	vector <unsigned char> readImage(bytesPerPixel * imageHeight* imageWidth);
	vector <unsigned char> readImageOriginal(bytesPerPixel * imageHeight* imageWidth);
	vector <unsigned char> readImageGray(imageHeight * imageWidth);
	vector <unsigned char> noise(bytesPerPixel * imageHeight * imageWidth);
	vector <unsigned char> denoise(bytesPerPixel * imageHeight * imageWidth);
	vector <unsigned char> histImage(bytesPerPixel * (imageHeight + (2 * N)) * (imageWidth + (2 * N)));


	double sigma = 2.0;

	if ((fopen_s(&imageRead, argv[1], "rb")))
	{
		cout << "Cannot open file" << argv[1] << endl;
		exit(1);
	}

	fread((unsigned char*)&readImage[0], sizeof(unsigned char), imageHeight* imageWidth*bytesPerPixel, imageRead);
	fclose(imageRead);


	if ((fopen_s(&imageReadOriginal, argv[2], "rb")))
	{
		cout << "Cannot open file" << argv[2] << endl;
		exit(1);
	}

	fread((unsigned char*)&readImageOriginal[0], sizeof(unsigned char), imageHeight* imageWidth*bytesPerPixel, imageRead);
	fclose(imageRead);


	unsigned char grayVal = 0;
	for (row_it = 0; row_it < imageHeight; row_it++)
	{
		for (col_it = 0; col_it < imageWidth; col_it++)
		{
			valFrom = row_it * imageWidth + col_it;
			valTo = valFrom;
			grayVal = (readImage[valFrom + 0] + readImage[valFrom + 1] + readImage[valFrom + 2]) / 3;
			readImageGray[valTo] = (unsigned char)grayVal;

		}
	}

	// Define the gaussian kernel 
	int kerRow_it = 0, kerCol_it = 0;
	double rtTwoPiSigmaSq = 1 / (sqrt(2 * M_PI) * sigma);
	double raiseToE;
	double weightAtPosR = 0, weightAtPosG = 0, weightAtPosB = 0;
	double weightAtPos;
	double sumWeight = 0.0;
	double twoSigmaSq = 2 * sigma * sigma;
	vector <double> gaussianKernel((2 * N + 1) *(2 * N + 1));
	int kernelSize = N;



	for (kerRow_it = 0; kerRow_it < (2 * N + 1); ++kerRow_it)
	{
		for (kerCol_it = 0; kerCol_it < (2 * N + 1); ++kerCol_it)
		{
			raiseToE = (((kerRow_it - N) * (kerRow_it - N)) + ((kerCol_it - N) * (kerCol_it - N))) / twoSigmaSq;
			weightAtPos = exp(-raiseToE) * rtTwoPiSigmaSq;
			gaussianKernel[kerRow_it * N + kerCol_it] = weightAtPos;
			sumWeight += weightAtPos;
		}
	}

	for (kerRow_it = 0; kerRow_it < (2 * N + 1); ++kerRow_it)
	{
		for (kerCol_it = 0; kerCol_it < (2 * N + 1); ++kerCol_it)
		{
			gaussianKernel[kerRow_it * N + kerCol_it] /= sumWeight;
		}
	}

	// Pad array 
	int pad_size = N;
	int minus_i, minus_j;


	int padded_nrows = imageHeight + 2 * pad_size;
	int padded_ncols = imageWidth + 2 * pad_size;

	for (int i = 0; i < padded_nrows; ++i) {
		for (int j = 0; j < padded_ncols; ++j) {
			for (int k = 0; k < bytesPerPixel; k++) {

				//reverse of i and j
				minus_i = padded_nrows - 1 - i;
				minus_j = padded_ncols - 1 - j;

				//cout << i << endl;

				//top edge
				if (i < pad_size) {
					if (j < pad_size) {
						histImage[(i * padded_ncols + j) * bytesPerPixel + k] = readImage[((pad_size - i) * imageHeight + (pad_size - j)) * bytesPerPixel + k];
					}
					else if (j >(imageWidth - 1 + pad_size)) {
						histImage[(i * padded_ncols + j) * bytesPerPixel + k] = readImage[((pad_size - i) * imageWidth + (j - pad_size - 2 * (pad_size - minus_j))) * bytesPerPixel + k];

					}
					else {
						histImage[(i * padded_ncols + j) * bytesPerPixel + k] = readImage[((pad_size - i) * imageWidth + (j - pad_size)) * bytesPerPixel + k];
					}
				}
				//bottom edge
				else if (i > imageHeight - 1 + pad_size) {
					if (j < pad_size) {
						histImage[(i * padded_ncols + j) * bytesPerPixel + k] = readImage[((i - pad_size - 2 * (pad_size - minus_i)) * imageWidth + (pad_size - j)) * bytesPerPixel + k];
					}
					else if (j > imageWidth + pad_size) {
						histImage[(i * padded_ncols + j) * bytesPerPixel + k] = readImage[((i - pad_size - 2 * (pad_size - minus_i)) * imageWidth + (j - pad_size - 2 * (pad_size - minus_j))) * bytesPerPixel + k];
					}
					else {
						histImage[(i * padded_ncols + j) * bytesPerPixel + k] = readImage[((i - pad_size - 2 * (pad_size - minus_i)) * imageWidth + (j - pad_size)) * bytesPerPixel + k];
					}
				}
				//left side
				else if (j < pad_size) {
					histImage[(i * padded_ncols + j) * bytesPerPixel + k] = readImage[((i - pad_size) * imageWidth + (pad_size - j)) * bytesPerPixel + k];
				}
				//right side
				else if (j > imageWidth - 1 + pad_size) {
					histImage[(i * padded_ncols + j) * bytesPerPixel + k] = readImage[((i - pad_size) * imageWidth + (j - pad_size - 2 * (pad_size - minus_j)))* bytesPerPixel + k];
				}
				//copy rest of array
				else {
					histImage[(i * padded_ncols + j) * bytesPerPixel + k] = readImage[((i - pad_size) * imageWidth + (j - pad_size)) * bytesPerPixel + k];
				}
			}
		}
	}
	vector<double> YRi((2 * N + 1)* (2 * N + 1));
	vector<double> YRj((2 * N + 1)* (2 * N + 1));
	vector<double> YGi((2 * N + 1)* (2 * N + 1));
	vector<double> YGj((2 * N + 1)* (2 * N + 1));
	vector<double> YBi((2 * N + 1)* (2 * N + 1));
	vector<double> YBj((2 * N + 1)* (2 * N + 1));

	double sumWeightR, sumWeightG, sumWeightB;
	weightAtPosR = 0;
	double test;
	int  neighbour = imageHeight + 2 * N;
	double diffrenceR, diffrenceG, diffrenceB, wmaxR, wmaxG, wmaxB, averageR, averageG, averageB, sweightR, sweightG, sweightB;
	int runTill = 2 * N + 1;
	int row1, col1, rowmin, rowmax, colmin, colmax;
	int count = 0;
	for (row_it = 0; row_it < imageHeight; ++row_it)
	{
		row1 = int(row_it + N);
		for (col_it = 0; col_it < imageWidth; ++col_it)
		{
			col1 = int(col_it + N);

			for (int i = 0; i < runTill; ++i)
			{
				for (int j = 0; j < runTill; ++j)
				{

					YRi[(i * runTill + j)] = histImage[((row1 - N + i) * imageWidth + (col1 - N + j)) * bytesPerPixel + 0];
					YGi[(i * runTill + j)] = histImage[((row1 - N + i) * imageWidth + (col1 - N + j)) * bytesPerPixel + 1];
					YBi[(i * runTill + j)] = histImage[((row1 - N + i) * imageWidth + (col1 - N + j)) * bytesPerPixel + 2];
				}
			}
			wmaxR = 0;
			wmaxG = 0;
			wmaxB = 0;
			averageR = 0;
			averageG = 0;
			averageB = 0;
			sumWeightR = 0;
			sumWeightG = 0;
			sumWeightB = 0;

			rowmin = max((row1 - N), N);
			rowmax = min((row1 + N), (imageHeight + N - 1));
			colmin = max((col1 - N), N);
			colmax = min((col1 + N), (imageWidth + N - 1));



			for (int row = rowmin; row < rowmax + 1; ++row)
			{
				for (int col = colmin; col < colmax + 1; ++col)
				{

					diffrenceR = 0;
					diffrenceG = 0;
					diffrenceB = 0;
					if (row == row1 && col == col1) {}
					else {

						for (int i = 0; i < runTill; ++i)
						{
							for (int j = 0; runTill < 5; ++j)
							{

								YRj[(i * runTill + j)] = histImage[((row + i - N) * imageWidth + (col - N + j)) * bytesPerPixel + 0];
								YGj[(i * runTill + j)] = histImage[((row + i - N) * imageWidth + (col - N + j)) * bytesPerPixel + 1];
								YBj[(i * runTill + j)] = histImage[((row + i - N) * imageWidth + (col - N + j)) * bytesPerPixel + 2];
								diffrenceR += gaussianKernel[(i * runTill + j)] * (YRi[(i * runTill + j)] - YRj[(i * runTill + j)]) *  (YRi[(i * runTill + j)] - YRj[(i * runTill + j)]);
								diffrenceG += gaussianKernel[(i * runTill + j)] * (YGi[(i * runTill + j)] - YGj[(i * runTill + j)]) *  (YGi[(i * runTill + j)] - YGj[(i * runTill + j)]);
								diffrenceB += gaussianKernel[(i * runTill + j)] * (YBi[(i * runTill + j)] - YBj[(i * runTill + j)]) *  (YBi[(i * runTill + j)] - YBj[(i * runTill + j)]);								//out << count++ << endl;

							}
						}
						weightAtPosR = exp(-(diffrenceR) / 4);
						weightAtPosG = exp(-(diffrenceG) / 4);
						weightAtPosB = exp(-(diffrenceB) / 4);

						if (weightAtPosR > wmaxR)
							wmaxR = weightAtPosR;
						sumWeightR += weightAtPosR;
						averageR += weightAtPosR * histImage[(row * neighbour + col) * bytesPerPixel + 0];


						if (weightAtPosG > wmaxG) wmaxG = weightAtPosG;
						sumWeightG += weightAtPosG;
						averageG += weightAtPosG * histImage[(row * neighbour + col) * bytesPerPixel + 1];

						if (weightAtPosB > wmaxB) wmaxB = weightAtPosB;
						sumWeightB += weightAtPosB;
						averageB += weightAtPosB * histImage[(row * neighbour + col) * bytesPerPixel + 2];



					}

					averageR += wmaxR * histImage[(row1 *neighbour + col1) * bytesPerPixel + 0];
					sumWeightR += wmaxR;

					averageG += wmaxG * histImage[(row1 * neighbour + col1) * bytesPerPixel + 1];
					sumWeightG += wmaxG;

					averageB += wmaxB * histImage[(row1 *neighbour + col1) * bytesPerPixel + 2];
					sumWeightB += wmaxB;


					if (sumWeightR > 0)

						denoise[(row_it * 512 + col_it) * bytesPerPixel + 0] = averageR / sumWeightR;

					else
						denoise[(row_it * 512 + col_it) * bytesPerPixel + 0] = readImage[(row_it * 512 + col_it) * bytesPerPixel + 0];

					if (sumWeightG > 0)

						denoise[(row_it * 512 + col_it) * bytesPerPixel + 1] = averageG / sumWeightG;

					else
						denoise[(row_it * 512 + col_it) * bytesPerPixel + 1] = readImage[(row_it * 512 + col_it) * bytesPerPixel + 1];
					if (sumWeightB > 0)

						denoise[(row_it * 512 + col_it) * bytesPerPixel + 2] = averageB / sumWeightB;

					else
						denoise[(row_it * 512 + col_it) * bytesPerPixel + 2] = readImage[(row_it * 512 + col_it) * bytesPerPixel + 2];

					//if (denoise[(row_it * 512 + col_it) * bytesPerPixel + 0] < 0)  denoise[(row_it * 512 + col_it) * bytesPerPixel + 0] = -9999;
					//noise[(row_it * 512 + col_it) * bytesPerPixel + pixel_it] = readImage[(row_it * 512 + col_it) * bytesPerPixel + pixel_it] - denoise[(row_it * 512 + col_it) * bytesPerPixel + pixel_it];
				}
			}

		}
	}


	//Calculate MSE
	double MSE[3] = { 0,0,0 };
	double PSNR[3] = { 0,0,0 };
	int newValue, oriValue;


	for (int k = 0; k < 3; k++)
	{
		for (int i = 0; i < imageHeight; i++)
		{
			for (int j = 0; j < imageWidth; j++)
			{

				newValue = denoise[(i * imageWidth + j) * 3 + k];
				oriValue = readImageOriginal[(i * imageWidth + j) * 3 + k];

				MSE[k] += (double)(newValue - oriValue) * (newValue - oriValue) / (imageHeight * imageWidth);
			}
		}
	}
	for (int k = 0; k < 3; k++)
	{
		PSNR[k] = 10 * log10((double)(255 * 255) / MSE[k]);
		printf("MSE:%lf  PSNR:%lf\n", MSE[k], PSNR[k]);
		out << k << " " << MSE[k] << " " << PSNR[k] << endl;
	}
	printf("******\n");

	if ((fopen_s(&imageWrite, argv[3], "wb")))
	{
		cout << "Cannot open output file " << endl;
	}

	fwrite((unsigned char*)&denoise[0], sizeof(unsigned char), imageHeight * imageWidth * bytesPerPixel, imageWrite);
	fclose(imageWrite);

}


