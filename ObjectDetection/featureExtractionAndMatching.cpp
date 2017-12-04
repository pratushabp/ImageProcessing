#include "featureExtractionAndMatching.h"

#include "opencv2/ml/ml.hpp"
using namespace cv::ml;
void featureExtraction(int argc, char * argv[])
{
	Mat img, imgGray, outimg;
	img = imread("bus.jpg", CV_LOAD_IMAGE_COLOR);
	Ptr<Feature2D> f2d = xfeatures2d::SURF::create();
	vector<KeyPoint> keypoints;
	f2d->detect(img, keypoints);
	drawKeypoints(img, keypoints, outimg);
	imwrite("extract_bus_surf.jpg", outimg);
	namedWindow("display", WINDOW_AUTOSIZE);
	imshow("dislay", outimg);
	waitKey(0);
	destroyAllWindows();

	img = imread("jeep.jpg", CV_LOAD_IMAGE_COLOR);
	//Ptr<Feature2D> f2d = xfeatures2d::SIFT::create();
	//vector<KeyPoint> keypoints;
	//f2d->detect(img, keypoints);
	drawKeypoints(img, keypoints, outimg);
	imwrite("extract_jeep_surf.jpg", outimg);
	namedWindow("display", WINDOW_AUTOSIZE);
	imshow("dislay", outimg);
	waitKey(0);
	destroyAllWindows();
}

void imageMatch(int argc, char * argv[])
{
	Mat img1, img2, descriptors1, descriptors2, outimg;
	img1 = imread("rav4_1.jpg", CV_LOAD_IMAGE_COLOR);
	img2 = imread("bus.jpg", CV_LOAD_IMAGE_COLOR);
	Ptr<Feature2D> f2d = xfeatures2d::SIFT::create();
	vector<KeyPoint> keypoints1, keypoints2;
	f2d->detect(img1, keypoints1);
	f2d->detect(img2, keypoints2);

	f2d->compute(img1, keypoints1, descriptors1);
	f2d->compute(img2, keypoints2,  descriptors2);

	vector<DMatch> matches, goodMatches;
	BFMatcher matcher;

	matcher.match(descriptors1, descriptors2, matches);

	// Compute the good matches

	double minDist = 100;  
	double maxDist = 0; 
	double dist;
	for(int  i = 0; i < descriptors1.rows; i++)
	{
		dist = matches[i].distance;
		minDist = dist < minDist ? dist : minDist;
	}
	
	//Threshold for SURF = 0.09
	//Threshold for SIFT = 0.02

	minDist = minDist > 0.2 ? 2 * minDist : 0.02;
	for (int i = 0; i < descriptors1.rows; i++)
	{
		if (matches[i].distance <= minDist)
		{
			goodMatches.push_back(matches[i]);
		}
	}
	
	drawMatches(img1, keypoints1, img2, keypoints2, goodMatches, outimg, Scalar::all(-1),Scalar::all(-1), vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);

	
	imwrite("matches_BF_sift_RavBus.jpg", outimg);
	namedWindow("display", WINDOW_AUTOSIZE);
	imshow("dislay", outimg); 
	waitKey(0);
	destroyAllWindows();

}

void bagOfWords(int argc, char * argv[])
{
	Ptr<Feature2D> f2d = xfeatures2d::SIFT::create();
	Ptr<Feature2D> SIFTDe = xfeatures2d::SiftDescriptorExtractor::create();
	//FlannBasedMatcher matcher;


	Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create("FlannBased");

	vector<KeyPoint> keyPointJeep, keyPointBus, keyPointRav1, keyPointRav2;
	Mat imgJeep, imgBus, imgRav1, imgRav2, fJeep, fBus, fRav1, fRav2, fbowJeep, fbowBus, fbowRav1, fbowRav2;
	//Train Imsges
	imgJeep = imread("Jeep.jpg", CV_8UC3);
	imgBus = imread("Bus.jpg", CV_8UC3);
	imgRav1 = imread("Rav4_1.jpg", CV_8UC3);
	imgRav2 = imread("Rav4_2.jpg", CV_8UC3);

	f2d->detect(imgJeep, keyPointJeep);
	f2d->detect(imgBus, keyPointBus);
	f2d->detect(imgRav1, keyPointRav1);
	f2d->detect(imgRav2, keyPointRav2);

	f2d->compute(imgJeep, keyPointJeep, fJeep);
	f2d->compute(imgBus, keyPointBus, fBus);
	f2d->compute(imgRav1, keyPointRav1, fRav1);
	f2d->compute(imgRav2, keyPointRav2, fRav2);
	fJeep.convertTo(fJeep, CV_32F);
	fBus.convertTo(fBus, CV_32F);
	fRav1.convertTo(fRav1, CV_32F);
	fRav2.convertTo(fRav2, CV_32F);
	int size = 8;
	
	cv::BOWKMeansTrainer bowTrainer(size, TermCriteria(TermCriteria::MAX_ITER + TermCriteria::EPS, 10, 0.001), 1, KMEANS_PP_CENTERS);
	
	bowTrainer.add(fJeep);
	bowTrainer.add(fBus);
	bowTrainer.add(fRav1);
	bowTrainer.add(fRav2);
	Mat codeBook = bowTrainer.cluster();
	vector<double> cb;
	
	ofstream out("codeWords.txt");
	
	Mat histResponse;
	cv::BOWImgDescriptorExtractor DE(SIFTDe, matcher );
	DE.setVocabulary(codeBook);
	DE.compute(imgJeep, keyPointJeep, fbowJeep);
	DE.compute(imgBus, keyPointBus, fbowBus);
	DE.compute(imgRav1, keyPointRav1, fbowRav1);
	DE.compute(imgRav2, keyPointRav2, fbowRav2);

	cout <<  codeBook << endl << "Codebook with 8 bins" <<endl;
	out << codeBook << endl << "Codebook with 8 bins" << endl;  
	//Train SVM

	int labels[4] = { 0,1,2,3 };
	Mat labelMat(4, 1, CV_32S, labels);

	Mat samples;
	samples.push_back(fbowJeep);
	samples.push_back(fbowBus);
	samples.push_back(fbowRav1);
	samples.push_back(fbowRav2);


		cout << "Train SVM.." << endl;
		Mat actualSamples;
		samples.convertTo(samples, CV_32F);
		

		Ptr<SVM> svm = SVM::create();
		svm->setType(SVM::C_SVC);
		svm->setKernel(SVM::LINEAR);
		svm->setTermCriteria(TermCriteria(TermCriteria::MAX_ITER, 100, 1e-6));
		svm->train(samples, ROW_SAMPLE, labelMat);

		Mat testImg = imread("Rav4_2.jpg", CV_8UC3);
		Mat histTest;
		vector<KeyPoint> keyTest;
		f2d->detect(testImg, keyTest);
		DE.compute(testImg, keyTest, histTest);

		int value = svm->predict(histTest);
		
		if (value == 0) cout << "Classified as label 0 - Jeep" << endl;
		if (value == 1) cout << "Classified as label 1 - Bus" << endl;
		if (value == 2) cout << "Classified as label 2 - RAV4_1" << endl;
		if (value == 3) cout << "Classified as label 3 - RAV4_2" << endl;

}

	

