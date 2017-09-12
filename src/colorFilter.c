/* ===============================================================
 * Copyright (C) 2017, Yirami  
 * Released under BSD License
 * All rights reserved.
 * 
 * Filename: colorFilter.c (compiled into mex)
 * F-Description: filter the true color image by color feature
 *
 *		input 	a MxNx3	array 	"image"
 * 			and 
 *		outputs a MxN 	matrix 	"filterImage"
 *
 * 		The calling syntax is:
 *			[filterImage] = colorFilter(image)
 * 
 * Author：Yirami
 * A-Description:
 * 		mailto:	  1446675043@qq.com
 * 		website:  https://yirami.github.io./
 * 
 * Date: 01-06-2017
 * Version: 1.0
 *  
 * New  Release:
 * N-Description:
 *
 * =============================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"

#define pi 3.14159265358979

/* The computational routine */
void colorFilter(unsigned char *image, size_t rows, size_t cols, unsigned char Threshold3, unsigned char Threshold2, unsigned char *filterImage)
{
// --> 
	size_t rBias = 0;
	size_t gBias = rows*cols-1;
	size_t bBias = 2*rows*cols-1;
	// "j" means cols-1
    for (size_t j=0; j<cols; j++) 
	{
		// and "i" means rows-1
		for (size_t i=0; i<rows; i++)
		{
			int thisR = (int)image[rBias+j*rows+i];
			int thisG = (int)image[gBias+j*rows+i];
			int thisB = (int)image[bBias+j*rows+i];
			if (thisR>=Threshold3 && thisG>=Threshold3 && thisB>=Threshold3)
			{
				/*if (abs(thisR-thisG)<(int)Threshold2 && abs(thisR-thisB)<(int)Threshold2 && abs(thisG-thisB)<(int)Threshold2)
				{*/
					filterImage[j*rows+i] = 255;
				/*}*/
			}
		}
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
// --> Parameters check

    // check for proper number of arguments
    if(nrhs!=1) 
	{
        mexErrMsgIdAndTxt("MyToolbox:colorFilter:nrhs","One input required.");
    }
    if(nlhs!=1) 
	{
        mexErrMsgIdAndTxt("MyToolbox:colorFilter:nlhs","One output required.");
    }
    // check the type of input argument(s)
    if( !mxIsUint8(prhs[0]) ) 
	{
        mexErrMsgIdAndTxt("MyToolbox:colorFilter:notUint8","Input matrix must be type Uint8.");
    }
    // check the shape of input argument(s)
    if(mxGetNumberOfDimensions(prhs[0])!=3) 
	{
        mexErrMsgIdAndTxt("MyToolbox:colorFilter:not3D","Input must be 3-D.");
    }
// --> Accept input parameters

    // create a pointer to the real data in the input matrix
    unsigned char *image = mxGetData(prhs[0]);
    // and get dimensions of the input matrix
	const mwSize *objSize = mxGetDimensions(prhs[0]);
	mwSize rows = objSize[0];
    mwSize cols = objSize[1];
	mwSize pages = objSize[2];
    // check the input image is true color
    if(pages!=3) 
	{
        mexErrMsgIdAndTxt("MyToolbox:colorFilter:notTrueColor","Input image must be true color.");
    }
// --> Parameters can be added to function interface

	// default parameters
	unsigned char Threshold3 = 150;
	unsigned char Threshold2 = 10;
	// check them
	
// --> Prepare data

	// creat the output matrix "filterImage" and pass the pointer
    plhs[0] = mxCreateNumericMatrix (rows,cols,mxUINT8_CLASS,mxREAL);
    unsigned char *filterImage = mxGetData(plhs[0]);
// --> Main computation

    // call the computational routine
    colorFilter(image,rows,cols,Threshold3,Threshold2,filterImage);
// --> Rest

	// release memory
	
}
