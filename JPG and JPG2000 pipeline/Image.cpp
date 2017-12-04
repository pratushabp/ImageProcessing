//*****************************************************************************
//
// Image.cpp : Defines the class operations on images
//
//
// Student implementation
//*****************************************************************************

#include "Image.h"
#include "math.h"
#define  M_PI 3.14159
#define Clamp255(n) ((n) > 255 ? 255 : ((n) < 0? 0 : (n)))


// Constructor and Desctructors
MyImage::MyImage() 
{
	Data = NULL;
	DCT = NULL;
	IDCT = NULL;
	DWT = NULL;
	IDWT = NULL;
	Width = -1;
	Height = -1;
	ImagePath[0] = 0;
}

MyImage::~MyImage()
{
	if ( Data )
		delete Data;
	if ( DCT )
		delete DCT;
	if ( IDCT )
		delete IDCT;
	if ( DWT )
		delete DWT;
	if ( IDWT )
		delete IDWT;

}

// Copy constructor
MyImage::MyImage( MyImage *otherImage)
{
	Height = otherImage->Height;
	Width  = otherImage->Width;
	Data   = new unsigned char[Width*Height*3];
	IDCT   = new unsigned char[Width*Height*3];
    DCT    = new float[Width*Height*3];
	IDWT   = new unsigned char[Width*Height*3];
    DWT    = new float[Width*Height*3];
	strcpy(otherImage->ImagePath, ImagePath );

	for ( int i=0; i<(Height*Width*3); i++ )
	{
		Data[i]	= otherImage->Data[i];
	}
}


/*char clamp(float t)
{
	t = floor(t+0.5);

	if(t > 
	return false;
}*/
// = operator overload
MyImage & MyImage::operator= (const MyImage &otherImage)
{
	Height = otherImage.Height;
	Width  = otherImage.Width;
	Data   = new unsigned char[Width*Height*3];           //3 bytes per pixel 
	IDCT   = new unsigned char[Width*Height*3]; 
	DCT    = new float[Width*Height*3];
	IDWT   = new unsigned char[Width*Height*3];
    DWT    = new float[Width*Height*3];
	strcpy( (char *)otherImage.ImagePath, ImagePath );

	for ( int i=0; i<(Height*Width*3); i++ )
	{
		Data[i]	= otherImage.Data[i];
	}
	
	return *this;

}


// MyImage::ReadImage
// Function to read the image given a path
bool MyImage::ReadImage()
{

	// Verify ImagePath
	if (ImagePath[0] == 0 || Width < 0 || Height < 0 )
	{
		fprintf(stderr, "Image or Image properties not defined");
		fprintf(stderr, "Usage is `Image.exe Imagefile w h`");
		return false;
	}
	
	// Create a valid output file pointer
	FILE *IN_FILE;
	IN_FILE = fopen(ImagePath, "rb");
	if ( IN_FILE == NULL ) 
	{
		fprintf(stderr, "Error Opening File for Reading");
		return false;
	}

	// Create and populate RGB buffers
	int i;
	unsigned char *Rbuf = new unsigned char[Height*Width]; 
	unsigned char *Gbuf = new unsigned char[Height*Width]; 
	unsigned char *Bbuf = new unsigned char[Height*Width]; 

	for (i = 0; i < Width*Height; i ++)
	{
		Rbuf[i] = fgetc(IN_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		Gbuf[i] = fgetc(IN_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		Bbuf[i] = fgetc(IN_FILE);
	}
	
	// Allocate Data structure and copy
	Data = new unsigned char[Width*Height*3];
	for (i = 0; i < Height*Width; i++)
	{
		Data[3*i]	= Bbuf[i];
		Data[3*i+1]	= Gbuf[i];
		Data[3*i+2]	= Rbuf[i];
	}

	// Clean up and return
	delete Rbuf;
	delete Gbuf;
	delete Bbuf;
	fclose(IN_FILE);

	return true;

}



// MyImage functions defined here
bool MyImage::WriteImage()
{
	// Verify ImagePath
	// Verify ImagePath
	if (ImagePath[0] == 0 || Width < 0 || Height < 0 )
	{
		fprintf(stderr, "Image or Image properties not defined");
		return false;
	}
	
	// Create a valid output file pointer
	FILE *OUT_FILE;
	OUT_FILE = fopen(ImagePath, "wb");
	if ( OUT_FILE == NULL ) 
	{
		fprintf(stderr, "Error Opening File for Writing");
		return false;
	}

	// Create and populate RGB buffers
	int i;
    unsigned char *Rbuf = new unsigned char[Height*Width]; 
	unsigned char *Gbuf = new unsigned char[Height*Width]; 
	unsigned char *Bbuf = new unsigned char[Height*Width]; 

	for (i = 0; i < Height*Width; i++)
	{
		Bbuf[i] = Data[3*i];
		Gbuf[i] = Data[3*i+1];
		Rbuf[i] = Data[3*i+2];
	}

	
	// Write data to file
	for (i = 0; i < Width*Height; i ++)
	{
		fputc(Rbuf[i], OUT_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		fputc(Gbuf[i], OUT_FILE);
	}
	for (i = 0; i < Width*Height; i ++)
	{
		fputc(Bbuf[i], OUT_FILE);
	}
	
	// Clean up and return
	delete Rbuf;
	delete Gbuf;
	delete Bbuf;
	fclose(OUT_FILE);

	return true;

}





//---------------DCT functions------------------
bool MyImage::RDCT()
{
    unsigned char* Rbuf = new unsigned char[Height * Width];
	float* Rbuf_DCT = new float[Height * Width];
	unsigned char** Rbuf_2d = NULL;
	Rbuf_2d = new unsigned char*[Height];
	float** Rbuf_2d_DCT = NULL;
	Rbuf_2d_DCT = 	new float*[Height];
	int CompCoeff = 0;
	//CompCoeff = floor((CoeffNum/4096) + 0.5);

	for( int i = 0; i < Height*Width; i++)
	{
		Rbuf[i] = Data[3 * i + 2 ]; 
	}

    for(int j = 0; j < Height; j++)
	{
		Rbuf_2d[j]  = new unsigned char[Width];
		Rbuf_2d_DCT[j] = new float[Width];

	}

	int offset = 0;
	for(int j = 0; j < Height; j++)
	{
		offset = j * Height;
      for( int i = 0; i < Width; i++)
	  {
		  Rbuf_2d[j][i] = Rbuf[offset + i];
	  }
	}
   for(int r = 0; r < 64; r++)
   {
	 for (int c = 0; c < 64; c++)
	 {
	   for (int v = 0; v < 8; v++)
	   {
		for( int u = 0; u < 8; u++)
		{
			//double sum = 0;
			Rbuf_2d_DCT[r * 8 + v][c * 8 + u] = 0.0;
			      double Cu = 1;
				  double Cv = 1;
				  if( u == 0)
				  {
					  Cu = sqrt(1/2.0);
				  }
				  if( v == 0)
				  {
					  Cv = sqrt(1/2.0); 
				  }  


			for(int j = 0; j < 8; j++)
			{
				for(int i = 0; i < 8; i++)
				{
				  
				  Rbuf_2d_DCT[r * 8 + v][c * 8 + u] +=  Rbuf_2d[r * 8 + j][c * 8 + i] * cos(((2*i + 1) * u * M_PI)/16)
					                                   * cos(((2*j + 1) * v * M_PI)/16);
                  
				}
			}
		
			Rbuf_2d_DCT[r * 8 + v][c * 8 + u]  *=  Cu * Cv * 0.25;
		
		}
      }
   }
 }

		offset = 0;
		for(int i = 0; i < Height; i++)
		{
			offset = Width * i;
			for( int j = 0; j < Width; j++)
			{
				Rbuf_DCT[offset + j] =  Rbuf_2d_DCT[i][j];
			}
		
		}
	
		for(int i = 0; i < Width*Height; i++)
	{
		DCT[3 * i + 2 ] = Rbuf_DCT[i];	
	}
        delete(Rbuf);
		delete(Rbuf_2d);
        delete(Rbuf_DCT);
		delete(Rbuf_2d_DCT);
		return false;
}








bool MyImage::GDCT()
{
    unsigned char* Gbuf = new unsigned char[Height * Width];
	float* Gbuf_DCT = new float[Height * Width];
	unsigned char** Gbuf_2d = NULL;
	Gbuf_2d = new unsigned char*[Height];
	float** Gbuf_2d_DCT = NULL;
	Gbuf_2d_DCT = 	new float*[Height];

	for( int i = 0; i < Height*Width; i++)
	{
		Gbuf[i] = Data[3 * i + 1 ]; 
	}

    for(int j = 0; j < Height; j++)
	{
	    Gbuf_2d[j]  = new unsigned char[Width];
		Gbuf_2d_DCT[j] = new float[Width];
	}

	int offset = 0;
	for(int j = 0; j < Height; j++)
	{
		offset = j * Height;
      for( int i = 0; i < Width; i++)
	  {
		  Gbuf_2d[j][i] = Gbuf[offset + i];
	  }
	}


   for(int r = 0; r < 64; r++)
	 for (int c = 0; c < 64; c++)
	   for (int v = 0; v < 8; v++)
		for( int u = 0; u < 8; u++)
		{
			Gbuf_2d_DCT[r * 8 + v][c * 8 + u] = 0.0;
			double Cu = 1;
				  double Cv = 1;
				  if( u == 0)
				  {
					  Cu = 1/sqrt(2.0);
				  }
				  if( v == 0)
				  {
					  Cv = 1/sqrt(2.0); 
				  }

			for(int j = 0; j < 8; j++)
			{
				for(int i = 0; i < 8; i++)
				{
				  
				  Gbuf_2d_DCT[r * 8 + v][c * 8 + u] +=  Gbuf_2d[r * 8 + j][c * 8 + i] * cos(((2*i + 1) * u * M_PI)/16)
					                                   * cos(((2*j + 1) * v * M_PI)/16);
                  
				}
			}
			Gbuf_2d_DCT[r * 8 + v][c * 8 + u] *= Cu * Cv * 0.25;
		}
	
		offset = 0;
		for(int i = 0; i < Height; i++)
		{
			offset = Width * i;
		for(int j = 0; j < Width; j++)
		{
            Gbuf_DCT[offset + j] = Gbuf_2d_DCT[i][j];

		}
		}

		for(int i = 0; i < Width*Height; i++)
	{
    
		DCT[3*i + 1 ] =  Gbuf_DCT[i];
		
	}

		delete(Gbuf);
		delete(Gbuf_2d);
        delete(Gbuf_DCT);
		delete(Gbuf_2d_DCT);
	    return false;
}


bool MyImage::BDCT()
{
    unsigned char* Bbuf = new unsigned char[Height * Width];
	float* Bbuf_DCT = new float[Height * Width];
	unsigned char** Bbuf_2d = NULL;
	Bbuf_2d = new unsigned char*[Height];
	float** Bbuf_2d_DCT = NULL;
	Bbuf_2d_DCT = 	new float*[Height];

	for( int i = 0; i < Height*Width; i++)
	{
		Bbuf[i] = Data[3 * i ]; 
	}

    for(int j = 0; j < Height; j++)
	{
		Bbuf_2d[j]  = new unsigned char[Width];
		Bbuf_2d_DCT[j] = new float[Width];
	}

	int offset = 0;
	for(int j = 0; j < Height; j++)
	{
		offset = j * Height;
      for( int i = 0; i < Width; i++)
	  {
		  Bbuf_2d[j][i] = Bbuf[offset + i];
	  }
	}


   for(int r = 0; r < 64; r++)
	 for (int c = 0; c < 64; c++)
	 {
	   for (int v = 0; v < 8; v++)
		for( int u = 0; u < 8; u++)
		{
			Bbuf_2d_DCT[r * 8 + v][c * 8 + u] = 0.0;
			      double Cu = 1;
				  double Cv = 1;
				  if( u == 0)
				  {
					  Cu = sqrt(1/2.0);
				  }
				  if( v == 0)
				  {
					  Cv = sqrt(1/2.0); 
				  }

			for(int j = 0; j < 8; j++)
			{
				 
				for(int i = 0; i < 8; i++)
				{
				 
				  Bbuf_2d_DCT[r * 8 + v][c * 8 + u] +=  Bbuf_2d[r * 8 + j][c * 8 + i] * cos(((2*i + 1) * u * M_PI)/16)
					                                   * cos(((2*j + 1) * v * M_PI)/16) ;
                  
				}
			}
			Bbuf_2d_DCT[r * 8 + v][c * 8 + u] *= Cu * Cv * 0.25;
		}
	 }
		offset = 0;
		for(int i = 0; i < Height; i++)
		{
			offset = Width * i;
		for(int j = 0; j < Width; j++)
		{
            Bbuf_DCT[offset + j] = Bbuf_2d_DCT[i][j];

		}
		}
		for(int i = 0; i < Width*Height; i++)
	{
    
		DCT[3*i] = Bbuf_DCT[i];
		
	}
    
	    delete(Bbuf);
		delete(Bbuf_2d);
        delete(Bbuf_DCT);
		delete(Bbuf_2d_DCT);


	return false;
}





bool MyImage::DCT_ZigZagScan(float **Buf, int CompCoeff) 
{
 
	int m = 8;
	int n =8;
	int h = 0;
    int i= 0, j = 0, up=1;
    bool turned = false;
    int d[2][2] = { { 1, -1 }, { -1, 1 } };
    int corner[2][4] = { { 1, 0, 0, 1 }, { 0, 1, 1, 0 } };
    //int CompCoeff = 2;

    while (i < m && j < n) 
	{
       if( CompCoeff != 0)
	   {
		CompCoeff--;
	   }
	   else
	   {
	   	*(*(Buf + i) + j) = 0.0; 
	   }
        if (i == 0 || j == 0 || i == m - 1 || j == n - 1) 
		{
            if (!turned) 
			{
                int k = 2 * (up * (j / (n - 1)) | (1 - up) * (i / (m - 1)));
                i += corner[up][k];
                j += corner[up][k + 1];
                turned = true;
                up = 1 - up;
                continue;
            } else
                turned = false;
        }
        i += d[up][0];
        j += d[up][1];
		
    
		}
	
return false;
}
bool MyImage::DWT_ZigZagScan(float **Buf, int CompCoeff) 
{
 
	int m = 512;
	int n =512;
	int h = 0;
    int i= 0, j = 0, up=1;
    bool turned = false;
    int d[2][2] = { { 1, -1 }, { -1, 1 } };
    int corner[2][4] = { { 1, 0, 0, 1 }, { 0, 1, 1, 0 } };
    //int CompCoeff = 4096;
	int count = 0;

    while (i < m && j < n) 
	{
       if( count<CompCoeff )
	   {
		//	CompCoeff--;
		count++;
	   }
	   else
	   {
	   	*(*(Buf + i) + j) = 0.0; 
		//count++;
	   }
        if (i == 0 || j == 0 || i == m - 1 || j == n - 1) 
		{
            if (!turned) 
			{
                int k = 2 * (up * (j / (n - 1)) | (1 - up) * (i / (m - 1)));
                i += corner[up][k];
                j += corner[up][k + 1];
                turned = true;
                up = 1 - up;
                continue;
            } else
                turned = false;
        }
        i += d[up][0];
        j += d[up][1];
		
    
		}
	
return false;
}
bool MyImage::RIDCT(int CompCoeff)
{
    float* Rbuf = new float[Height * Width];
	unsigned char* Rbuf_DCT = new unsigned char[Height * Width];
	float** Rbuf_2d = NULL;
	Rbuf_2d = new float*[Height];
	float** Rbuf_2d_Comp = NULL;
	Rbuf_2d_Comp = new float*[Height];
	float** Rbuf_2d_DCT = NULL;
	Rbuf_2d_DCT = 	new float*[Height];

	CompCoeff = floor((CompCoeff/4096)+0.5);
	

	for( int i = 0; i < Height*Width; i++)
	{
		Rbuf[i] = DCT[3 * i + 2 ]; 
	}

    for(int j = 0; j < Height; j++)
	{
		Rbuf_2d[j]  = new float[Width];
		Rbuf_2d_DCT[j] = new float[Width];
		Rbuf_2d_Comp[j] = new float[Width];
	}


	int offset = 0;
	for(int j = 0; j < Height; j++)
	{
		offset = j * Height;
      for( int i = 0; i < Width; i++)
	  {
		  Rbuf_2d[j][i] = Rbuf[offset + i];
	  }
	}

	for(int r = 0; r < 64; r++)
	{
	 for (int c = 0; c < 64; c++)
	 {
		 for( int i = 0; i < 8; i++)
		 {
			 for(int j = 0; j < 8; j++)
			 {

				 Rbuf_2d_Comp[i][j] = Rbuf_2d[r * 8 + i][c * 8 + j]; 
   			     
	         }
		 }
		 DCT_ZigZagScan(Rbuf_2d_Comp,CompCoeff);
		 for( int i = 0; i < 8; i++)
		 {
			 for(int j = 0; j < 8; j++)
			 {
                 Rbuf_2d[r * 8 + i][c * 8 + j] = Rbuf_2d_Comp[i][j];
			 }
		 }
	   }
	}



   for(int r = 0; r < 64; r++)
	 for (int c = 0; c < 64; c++)
	   for (int j = 0; j < 8; j++)
		for( int i = 0; i < 8; i++)
		{
			Rbuf_2d_DCT[r * 8 + j][c * 8 + i] = 0.0;
			     
			for(int v = 0; v < 8; v++)
			{
				for(int u = 0; u < 8; u++)
				{
				  double Cu = 1;
				  double Cv = 1;
				  if( u == 0)
				  {
					  Cu = sqrt(1/2.0);
				  }
				  if( v == 0)
				  {
					  Cv = sqrt(1/2.0); 
				  }  


				  Rbuf_2d_DCT[r * 8 + j][c * 8 + i] +=  Rbuf_2d[r * 8 + v][c * 8 + u] * cos(((2*i + 1) * u * M_PI)/16)
					                                   * cos(((2*j + 1) * v * M_PI)/16) * Cu * Cv ;
                  
				}
			}
			 Rbuf_2d_DCT[r * 8 + j][c * 8 + i] *= 0.25;
		}

		offset = 0;
		for(int i = 0; i < Height; i++)
		{
			offset = Width * i;
			for( int j = 0; j < Width; j++)
			{
				Rbuf_DCT[offset + j] = (unsigned char)(Clamp255((floor(Rbuf_2d_DCT[i][j]+ 0.5))));
			}
		
		}
	
		for(int i = 0; i < Width*Height; i++)
	{
    
		IDCT[3*i + 2 ] = Rbuf_DCT[i];
		
	}


		delete(Rbuf);
		delete(Rbuf_2d);
        delete(Rbuf_2d_Comp);
		delete(Rbuf_2d_DCT);

	return false;
}



bool MyImage::GIDCT(int CompCoeff)
{
    float* Gbuf = new float[Height * Width];
	unsigned char* Gbuf_DCT = new unsigned char[Height * Width];
	float** Gbuf_2d = NULL;
	Gbuf_2d = new float*[Height];
	float** Gbuf_2d_Comp = NULL;
	Gbuf_2d_Comp = new float*[Height];
	float** Gbuf_2d_DCT = NULL;
	Gbuf_2d_DCT = 	new float*[Height];


	CompCoeff = floor((CompCoeff/4096)+0.5);

	for( int i = 0; i < Height*Width; i++)
	{
		Gbuf[i] = DCT[3 * i + 1 ]; 
	}

    for(int j = 0; j < Height; j++)
	{
		Gbuf_2d[j]  = new float[Width];
		Gbuf_2d_Comp[j] = new float[Width];
		Gbuf_2d_DCT[j] = new float[Width];
	}

	int offset = 0;
	for(int j = 0; j < Height; j++)
	{
		offset = j * Height;
      for( int i = 0; i < Width; i++)
	  {
		  Gbuf_2d[j][i] = Gbuf[offset + i];
	  }
	}

   for(int r = 0; r < 64; r++)
	{
	 for (int c = 0; c < 64; c++)
	 {
		 for( int i = 0; i < 8; i++)
		 {
			 for(int j = 0; j < 8; j++)
			 {

				 Gbuf_2d_Comp[i][j] = Gbuf_2d[r * 8 + i][c * 8 + j]; 
   			     
	         }
		 }
		 DCT_ZigZagScan(Gbuf_2d_Comp, CompCoeff);
		 for( int i = 0; i < 8; i++)
		 {
			 for(int j = 0; j < 8; j++)
			 {
                 Gbuf_2d[r * 8 + i][c * 8 + j] = Gbuf_2d_Comp[i][j];
			 }
		 }
	   }
	}


   for(int r = 0; r < 64; r++)
	 for (int c = 0; c < 64; c++)
	   for (int j = 0; j < 8; j++)
		for( int i = 0; i < 8; i++)
		{
			Gbuf_2d_DCT[r * 8 + j][c * 8 + i] = 0.0;
			     
			for(int v = 0; v < 8; v++)
			{
				for(int u = 0; u < 8; u++)
				{
				  double Cu = 1;
				  double Cv = 1;
				  if( u == 0)
				  {
					  Cu = sqrt(1/2.0);
				  }
				  if( v == 0)
				  {
					  Cv = sqrt(1/2.0); 
				  }  


				  Gbuf_2d_DCT[r * 8 + j][c * 8 + i] +=  Gbuf_2d[r * 8 + v][c * 8 + u] * cos(((2*i + 1) * u * M_PI)/16)
					                                   * cos(((2*j + 1) * v * M_PI)/16) * Cu * Cv ;
                  
				}
			}
			 Gbuf_2d_DCT[r * 8 + j][c * 8 + i] *= 0.25;
		}

		offset = 0;
		for(int i = 0; i < Height; i++)
		{
			offset = Width * i;
			for( int j = 0; j < Width; j++)
			{
				Gbuf_DCT[offset + j] = (unsigned char)(Clamp255((floor(Gbuf_2d_DCT[i][j]+ 0.5))));
			}
		
		}
	
		for(int i = 0; i < Width*Height; i++)
	{
    
		IDCT[3*i + 1 ] = Gbuf_DCT[i];
		
	}



   delete(Gbuf);
		delete(Gbuf_2d);
        delete(Gbuf_2d_Comp);
		delete(Gbuf_2d_DCT);
	return false;
}







bool MyImage::BIDCT(int CompCoeff)
{
    float* Bbuf = new float[Height * Width];
	unsigned char* Bbuf_DCT = new unsigned char[Height * Width];
	float** Bbuf_2d = NULL;
	Bbuf_2d = new float*[Height];
	float** Bbuf_2d_Comp = NULL;
	Bbuf_2d_Comp = new float*[Height];
	float** Bbuf_2d_DCT = NULL;
	Bbuf_2d_DCT = 	new float*[Height];


    CompCoeff = floor((CompCoeff/4096)+0.5);

	for( int i = 0; i < Height*Width; i++)
	{
		Bbuf[i] = DCT[3 * i ]; 
	}

    for(int j = 0; j < Height; j++)
	{
		Bbuf_2d[j]  = new float[Width];
		Bbuf_2d_DCT[j] = new float[Width];
		Bbuf_2d_Comp[j] = new float[Width];
	}

	int offset = 0;
	for(int j = 0; j < Height; j++)
	{
		offset = j * Height;
      for( int i = 0; i < Width; i++)
	  {
		  Bbuf_2d[j][i] = Bbuf[offset + i];
	  }
	}
 
	for(int r = 0; r < 64; r++)
	{
	 for (int c = 0; c < 64; c++)
	 {
		 for( int i = 0; i < 8; i++)
		 {
			 for(int j = 0; j < 8; j++)
			 {

				 Bbuf_2d_Comp[i][j] = Bbuf_2d[r * 8 + i][c * 8 + j]; 
   			     
	         }
		 }
		 DCT_ZigZagScan(Bbuf_2d_Comp, CompCoeff);
		 for( int i = 0; i < 8; i++)
		 {
			 for(int j = 0; j < 8; j++)
			 {
                 Bbuf_2d[r * 8 + i][c * 8 + j] = Bbuf_2d_Comp[i][j];
			 }
		 }
	   }
	}


   for(int r = 0; r < 64; r++)
	 for (int c = 0; c < 64; c++)
	   for (int j = 0; j < 8; j++)
		for( int i = 0; i < 8; i++)
		{
			Bbuf_2d_DCT[r * 8 + j][c * 8 + i] = 0.0;
			     
			for(int v = 0; v < 8; v++)
			{
				for(int u = 0; u < 8; u++)
				{
				  double Cu = 1;
				  double Cv = 1;
				  if( u == 0)
				  {
					  Cu = sqrt(1/2.0);
				  }
				  if( v == 0)
				  {
					  Cv = sqrt(1/2.0); 
				  }  


				  Bbuf_2d_DCT[r * 8 + j][c * 8 + i] +=  Bbuf_2d[r * 8 + v][c * 8 + u] * cos(((2*i + 1) * u * M_PI)/16)
					                                   * cos(((2*j + 1) * v * M_PI)/16) * Cu * Cv ;
                  
				}
			}
			 Bbuf_2d_DCT[r * 8 + j][c * 8 + i] *= 0.25;
		}

		offset = 0;
		for(int i = 0; i < Height; i++)
		{
			offset = Width * i;
			for( int j = 0; j < Width; j++)
			{
				Bbuf_DCT[offset + j] = (unsigned char)(Clamp255((floor(Bbuf_2d_DCT[i][j]+ 0.5))));
			}
		
		}
	
		for(int i = 0; i < Width*Height; i++)
	{
    
		IDCT[3*i ] = Bbuf_DCT[i];
		
	}


		delete(Bbuf);
		delete(Bbuf_2d);
        delete(Bbuf_2d_Comp);
		delete(Bbuf_2d_DCT);
	return false;
}











//-------------------DWT-------------------------------//

bool MyImage::RDWT()
{
    unsigned char* Rbuf = new unsigned char[Height * Width];
    float  ** Rbuf_2d = NULL;
	Rbuf_2d = new float*[Height];
	float* Rbuf_DWT = new float[Height * Width];

	int offset = 0;

	for( int i = 0; i < Height * Width; i++)
	{
		Rbuf[i] = Data[3*i+2];
	}

	for(int i = 0; i < Height; i++)
	{
		Rbuf_2d[i] = new float[Width];

	}

	for(int j = 0; j < Height; j++)
	{
		offset = Width*j;
		for(int i = 0; i < Width; i++)
		{
			Rbuf_2d[j][i] =(float) Rbuf[i + offset]; 
		}
	}

    DWT_SubBand(9,Rbuf_2d);

	for(int j = 0; j < Height; j++)
	{
		offset = Width * j;
		for(int i = 0; i < Width; i++)
		{

			Rbuf_DWT[offset + i] = Rbuf_2d[j][i];

		}
	}

	for(int i = 0; i < Height * Width; i++)
	{
		DWT[3 * i + 2] = Rbuf_DWT[i]; 
		//IDWT[3 * i + 2] =(unsigned char) Clamp255(floor(Rbuf_DWT[i]+0.5));
	}


	delete(Rbuf);
	delete(Rbuf_2d);
	delete(Rbuf_DWT);
	return false;
}








bool MyImage::RIDWT(int Coeff)
{
    float* Rbuf = new float[Height * Width];
    float  ** Rbuf_2d = NULL;
	Rbuf_2d = new float*[Height];
	float* Rbuf_DWT = new float[Height * Width];

	int offset = 0;

	for( int i = 0; i < Height * Width; i++)
	{
		Rbuf[i] = DWT[3*i + 2];
	}

	for(int i = 0; i < Height; i++)
	{
		Rbuf_2d[i] = new float[Width];

	}

	for(int j = 0; j < Height; j++)
	{
		offset = Width*j;
		for(int i = 0; i < Width; i++)
		{
			Rbuf_2d[j][i] = Rbuf[i + offset]; 
		}
	}
   
	 DWT_ZigZagScan(Rbuf_2d,Coeff);
     IDWT_SubBand(9,Rbuf_2d);
   
	for(int j = 0; j < Height; j++)
	{
		offset = Width * j;
		for(int i = 0; i < Width; i++)
		{

			Rbuf_DWT[offset + i] =Rbuf_2d[j][i];

		}
	}

	for(int i = 0; i < Height * Width; i++)
	{
		IDWT[3 * i + 2] = (unsigned char)(Clamp255((floor(Rbuf_DWT[i]+0.5)))); 
	}


	delete(Rbuf);
	delete(Rbuf_2d);
	delete(Rbuf_DWT);
	return false;
}







bool MyImage::GDWT()
{
    unsigned char* Gbuf = new unsigned char[Height * Width];
    float  ** Gbuf_2d = NULL;
	Gbuf_2d = new float*[Height];
	float* Gbuf_DWT = new float[Height * Width];

	int offset = 0;

	for( int i = 0; i < Height * Width; i++)
	{
		Gbuf[i] = Data[3 * i + 1];
	}

	for(int i = 0; i < Height; i++)
	{
		Gbuf_2d[i] = new float[Width];

	}

	for(int j = 0; j < Height; j++)
	{
		offset = Width*j;
		for(int i = 0; i < Width; i++)
		{
			Gbuf_2d[j][i] =(float) Gbuf[i + offset]; 
		}
	}

    DWT_SubBand(9,Gbuf_2d);

	for(int j = 0; j < Height; j++)
	{
		offset = Width * j;
		for(int i = 0; i < Width; i++)
		{

			Gbuf_DWT[offset + i] = Gbuf_2d[j][i];

		}
	}

	for(int i = 0; i < Height * Width; i++)
	{
		DWT[3 * i + 1] = Gbuf_DWT[i]; 
		//IDWT[3 * i + 1] = (unsigned char)Clamp255(Gbuf_DWT[i]);
	}

	delete(Gbuf);
	delete(Gbuf_2d);
	delete(Gbuf_DWT);

	return false;
}








bool MyImage::GIDWT(int Coeff)
{
    float* Gbuf = new float[Height * Width];
    float  ** Gbuf_2d = NULL;
	Gbuf_2d = new float*[Height];
	float* Gbuf_DWT = new float[Height * Width];

	int offset = 0;

	for( int i = 0; i < Height * Width; i++)
	{
		Gbuf[i] = DWT[3*i + 1];
	}

	for(int i = 0; i < Height; i++)
	{
		Gbuf_2d[i] = new float[Width];

	}

	for(int j = 0; j < Height; j++)
	{
		offset = Width*j;
		for(int i = 0; i < Width; i++)
		{
			Gbuf_2d[j][i] =(float) Gbuf[i + offset]; 
		}
	}

	DWT_ZigZagScan(Gbuf_2d,Coeff); 
    IDWT_SubBand(9,Gbuf_2d);
	

	for(int j = 0; j < Height; j++)
	{
		offset = Width * j;
		for(int i = 0; i < Width; i++)
		{

			Gbuf_DWT[offset + i] =Gbuf_2d[j][i];

		}
	}

	for(int i = 0; i < Height * Width; i++)
	{
		IDWT[3 * i + 1] = (unsigned char)(Clamp255((floor(Gbuf_DWT[i]+0.5)))); 
	}

	delete(Gbuf);
	delete(Gbuf_2d);
	delete(Gbuf_DWT);

	return false;
}







bool MyImage::BDWT()
{
    unsigned char* Bbuf = new unsigned char[Height * Width];
    float  ** Bbuf_2d = NULL;
	Bbuf_2d = new float*[Height];
	float* Bbuf_DWT = new float[Height * Width];

	int offset = 0;

	for( int i = 0; i < Height * Width; i++)
	{
		Bbuf[i] = Data[3*i];
	}

	for(int i = 0; i < Height; i++)
	{
		Bbuf_2d[i] = new float[Width];

	}

	for(int j = 0; j < Height; j++)
	{
		offset = Width*j;
		for(int i = 0; i < Width; i++)
		{
			Bbuf_2d[j][i] =(float) Bbuf[i + offset]; 
		}
	}

    DWT_SubBand(9,Bbuf_2d);

	for(int j = 0; j < Height; j++)
	{
		offset = Width * j;
		for(int i = 0; i < Width; i++)
		{

			Bbuf_DWT[offset + i] = Bbuf_2d[j][i];

		}
	}

	for(int i = 0; i < Height * Width; i++)
	{
		DWT[3 * i] = Bbuf_DWT[i]; 
		//IDWT[3 * i ] = (unsigned char)Clamp255(Bbuf_DWT[i]);
	}
	delete(Bbuf);
	delete(Bbuf_2d);
	delete(Bbuf_DWT);



	return false;
}








bool MyImage::BIDWT(int Coeff)
{
    float* Bbuf = new float[Height * Width];
    float  ** Bbuf_2d = NULL;
	Bbuf_2d = new float*[Height];
	float* Bbuf_DWT = new float[Height * Width];

	int offset = 0;

	for( int i = 0; i < Height * Width; i++)
	{
		Bbuf[i] = DWT[3*i];
	}

	for(int i = 0; i < Height; i++)
	{
		Bbuf_2d[i] = new float[Width];

	}

	for(int j = 0; j < Height; j++)
	{
		offset = Width*j;
		for(int i = 0; i < Width; i++)
		{
			Bbuf_2d[j][i] = Bbuf[i + offset]; 
		}
	}
	
	DWT_ZigZagScan(Bbuf_2d, Coeff);
    IDWT_SubBand(9,Bbuf_2d);
	 

	for(int j = 0; j < Height; j++)
	{
		offset = Width * j;
		for(int i = 0; i < Width; i++)
		{

			Bbuf_DWT[offset + i] = Bbuf_2d[j][i];

		}
	}

	for(int i = 0; i < Height * Width; i++)
	{
		IDWT[3 * i] = (unsigned char)(Clamp255((floor(Bbuf_DWT[i]+0.5)))); 
	}
 

	return false;
}


bool MyImage::DWT_SubBand(int n, float **Buf)
{
    
	
		DWT_RC(9,Buf);
	return false;
}

bool MyImage::DWT_RC(int n, float **Buf)
{
 

  //int coeffs;
  float*Tbuf = new float[Height * Width];
 // coeffs = pow(2,float(n));
  for (int i = 0;  i < 512;  i++)
	      {
			  for(int j=0;j<Height;j++)
			  {
				Tbuf[j]=Buf[i][j];
			  }
			DWT_1d_SubBand(n,Tbuf);
			for(int j=0;j<512;j++)
			  {
				Buf[i][j]=Tbuf[j];
			  }

    }

  //Tbuf = new float[coeffs];

  
  for (int j = 0;  j < Height;  j++)
    {
      for (int i = 0;  i < Width;  i++)
	  {
 	    Tbuf[i] = Buf[i][j];
 	  }
	  DWT_1d_SubBand(n,Tbuf);
      for (int i = 0;  i < Height;  i++)
	  {
		 Buf[i][j] = Tbuf[i];
	  }
    }




  //free (Tbuf);   
    delete Tbuf;
  return false;
}



bool MyImage::DWT_1d_SubBand(int n, float *Buf)
{
  

  for (int i = n-1;  i >= 0;  i--)
    {
      DWT_1d(i+1,Buf);
    }

	return false;
}

bool MyImage::DWT_1d(int n, float *Buf)

{
  int i;
  int coeffs;
  float * LP= NULL;
  float * HP = NULL;

  coeffs = pow(2,float(n));
  LP = new float[coeffs/2];
  HP= new float[coeffs/2];
  
  for (i = 0;  i < coeffs/2;  i++)
    {
      LP[i] = ((Buf[2*i] + Buf[2*i+1]) / 2.0);
      HP[i] = ((Buf[2*i] - Buf[2*i+1]) / 2.0);
    }

  for (i = 0;  i < coeffs/2;  i++)
    {
      Buf[i] = LP[i];
      Buf[i + coeffs/2] = HP[i];
    }

  free (LP);   
  LP = NULL;
  free (HP);   
  HP = NULL;
	return false;
}







bool MyImage::IDWT_SubBand(int n, float **Buf)
{
    //for( int i = 1; i <= n; i++)
	{
	  IDWT_RC(9,Buf);
	}
	return false;
}


bool MyImage::IDWT_RC(int n, float **Buf)
{
 

  int size = 512;
  
  //coeffs = pow(2,float(n));
  


  float*Tbuf = new float[Height * Width];

   for (int j = 0;  j < Height;  j++)
    {
      for (int i = 0;  i < Width;  i++)
	  Tbuf[i] = Buf[i][j];
	  {
		IDWT_1d_SubBand(n,Tbuf);
	  }
      for (int i = 0;  i < Height;  i++)
	  Buf[i][j] = Tbuf[i];
    }


  for (int i = 0;  i < Height;  i++)
    {
		for (int j = 0;  j < Width;  j++)
		{
			Tbuf[j]=Buf[i][j];
		}
      IDWT_1d_SubBand(n,Tbuf);
	  for (int j = 0;  j < Width;  j++)
		{
			Buf[i][j]=Tbuf[j];
		}
    }

  
  
  free (Tbuf);   
  Tbuf= NULL;
  return false;
}



bool  MyImage::IDWT_1d_SubBand(int n, float *Buf)
{
  

  for (int i = 1;  i <= n;  i++)
    {
      IDWT_1d(i,Buf);
    }

	return false;
}

bool MyImage::IDWT_1d(int n,float *Buf)
{
  int i;
  int coeffs;
  float *IDWT = NULL;

  coeffs = pow(2,float(n));
  IDWT = new float[coeffs];
  
  for (i = 0;  i < coeffs/2;  i++)
    {  
      IDWT[2*i]   = Buf[i] + Buf[i + coeffs/2];
      IDWT[2*i+1] = Buf[i] - Buf[i + coeffs/2]; 
    }

  for (i = 0;  i < coeffs;  i++)
    {
      Buf[i] = IDWT[i];
    }

  free (IDWT);   
  IDWT = NULL;
  return false;
}



//----------------------------------------------------//
bool MyImage::DCT_IDCT(int Coeff_t)
{
	 BDCT();
	 GDCT();
     RDCT();
	 BIDCT(Coeff_t);
	 GIDCT(Coeff_t);
	 RIDCT(Coeff_t);

	return false;
}
bool MyImage::DWT_IDWT(int Coeff)
{
	RDWT();
 	RIDWT(Coeff);
	GDWT();
	GIDWT(Coeff);
	BDWT();
	BIDWT(Coeff);
	return false;
}

