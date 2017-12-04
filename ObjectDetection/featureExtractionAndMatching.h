#pragma once


#include "myImage.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2\xfeatures2d\nonfree.hpp"
#include "opencv2\xfeatures2d.hpp"
#include <iostream>

void featureExtraction(int argc, char*argv[]);
void imageMatch(int argc, char*argv[]);
void bagOfWords(int argc, char*argv[]);
