#include <stdio.h>
#include <iostream>
#include <conio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
	cout << "Digital Image Processing" << endl;
	cout << "Assignment 1 " << endl;
	cout << "Pratusha Bhuvana Prasad" << endl;
	if (!(strcmp(argv[1], "crop")))
		crop(argc - 1, argv + 1);
	if (!(strcmp(argv[1], "resize")))
		resize(argc - 1, argv + 1);
	if (!(strcmp(argv[1], "cmyk")))
		cmyk(argc - 1, argv + 1);
	if (!(strcmp(argv[1], "hsl")))
		hsl(argc - 1, argv + 1);
	if (!(strcmp(argv[1], "histEqualisationGrayBucket")))
		histEqualisationGrayBucket(argc - 1, argv + 1);
	if (!(strcmp(argv[1], "histEqualisationGrayTF")))
		histEqualisationGrayTF(argc - 1, argv + 1);
	if (!(strcmp(argv[1], "histEqualisationColorBucket")))
		histEqualisationColorBucket(argc - 1, argv + 1);
	if (!(strcmp(argv[1], "histEqualisationColorTF")))
		histEqualisationColorTF(argc - 1, argv + 1);
	if (!(strcmp(argv[1], "specialEffect")))
		specialEffect(argc - 1, argv + 1);
	if (!(strcmp(argv[1], "histTransform")))
		histTransform(argc - 1, argv + 1);
	if (!(strcmp(argv[1], "mixNoise")))
		mixNoise(argc - 1, argv + 1);
	if (!(strcmp(argv[1], "nonLocalMean")))
		nonLocalMean(argc - 1, argv + 1);

	else
	{
		cerr << "Usage:" << endl;
		cerr << "crop" << endl;
		cerr << "resize" << endl;
		cerr << "cmyk" << endl;
		cerr << "hsl" << endl;
		cerr << "histEqualisationGray" << endl;
		cerr << "histEqualisationColor" << endl;
		cerr << "specialEffect" << endl;
		cerr << "histTransform" << endl;
		cerr << "mixNoise" << endl;
		cerr << "nonLocalMean" << endl;
		cerr << "bm3d" << endl;

	}
}