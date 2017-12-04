//*****************************************************************************
//
// Image.h : Defines the class operations on images
//
//*****************************************************************************

#ifndef IMAGE_DISPLAY
#define IMAGE_DISPLAY


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "resource.h"
#include "afxwin.h"

// C RunTime Header Files
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>


// Class structure of Image 
// Use to encapsulate an RGB image
class MyImage 
{

private:
	int		Width;					// Width of Image
	int		Height;					// Height of Image
	char	ImagePath[_MAX_PATH];	// Image location
	unsigned char* Data;					// RGB data of the image
	float* DCT;
	float* DCT_Compressed;
	unsigned char* IDCT;
	float* DWT;
	unsigned char* IDWT;
public:
	// Constructor
	MyImage();
	// Copy Constructor
	MyImage::MyImage( MyImage *otherImage);
	// Destructor
	~MyImage();

	// operator overload
	MyImage & operator= (const MyImage & otherImage);

	// Reader & Writer functions
	void	setWidth( const int w)  { Width = w; }; 
	void	setHeight(const int h) { Height = h; }; 
	void	setImageData( const char *img ) { Data = (unsigned char *)img; };
	void	setImagePath( const char *path) { strcpy(ImagePath, path); }
	int		getWidth() { return Width; };
	int		getHeight() { return Height; };
	unsigned char*	getImageData() { return Data; };
	unsigned char*	getDCTData() { return (IDCT); };
	
	unsigned char*	getDWTData() { return (IDWT); };
	char*	getImagePath() { return ImagePath; }

	// Input Output operations
	bool	ReadImage();
	bool	WriteImage();
    bool    DCT_ZigZagScan(float**, int);
	bool    DWT_ZigZagScan(float **Buf, int);
	// Modifications
//	bool	DCT();
	bool    BDCT();
	bool    GDCT();
	bool    RDCT();
	bool    RIDCT(int);
	bool    GIDCT(int);
	bool    BIDCT(int);
	bool    DWT_1d(int n, float*Buf);
    bool    DWT_1d_SubBand(int n, float*Buf);
	bool    DWT_RC(int n, float**Buf);
	bool    DWT_SubBand(int n, float**Buf);
	bool    IDWT_1d(int n, float*Buf);
    bool    IDWT_1d_SubBand(int n, float*Buf);
	bool    IDWT_RC(int n, float**Buf);
	bool    IDWT_SubBand(int n, float**Buf);
    bool    RDWT();
	bool    RIDWT(int);
    bool    GDWT();
	bool    GIDWT(int);
	bool    BDWT();
	bool    BIDWT(int);
	bool    DCT_IDCT(int);
    bool    DWT_IDWT(int);


};

#endif //IMAGE_DISPLAY
