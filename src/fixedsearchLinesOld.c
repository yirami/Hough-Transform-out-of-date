/* ===============================================================
 * Copyright (C) 2017, Yirami  
 * Released under BSD License
 * All rights reserved.
 * 
 * Filename: fixedsearchLines.c (compiled into mex)
 * F-Description: search lines for various Hough Transform
 *
 *		inputs 	a MxN matrix "img"
 *			 	a PxQ matrix "vote"
 *			 	a 	  vector "thetaSeq"
 *			 	a 	  vector "rhoSeq"
 *				a 	  scalar "linesMax"
 *				a 	  scalar "lineGap"
 *			 	a 	  scalar "angleS"
 *			 	a	  scalar "angleE"
 * 			and 
 *		output	a 6xR matrix "linesInfo"
 *			output format	-> 	start	point(r)
 *							-> 	start 	point(c)
 *							-> 	end 	point(r)
 *							-> 	end 	point(c)
 *							-> 	line 	theta
 *							-> 	line	rho
 *
 * 		The calling syntax is:
 *			[linesInfo] = fixedsearchLines(img,vote,thetaSeq,rhoSeq,
 *				... linesMax,lineGap,angleS,angleE)
 * 
 * Author：Yirami
 * A-Description:
 * 		mailto:	  1446675043@qq.com
 * 		website:  https://yirami.github.io./
 * 
 * Date: 26-05-2017
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

/* Procedure of searching peaks */
void peaks(unsigned int *vote, size_t m, size_t n, unsigned int linesMax, double *thetaSeq, double*rhoSeq, double thetaStep,double rhoStep, double *peaksInfo, size_t pre_subs_left, size_t pre_subs_right)
{
	// built-in parameters
	double suppressionVector[2] = {1,2};
	// prepare data
	size_t suppVector[2];
	suppVector[0] = (size_t)ceil(thetaStep*suppressionVector[0]);
	suppVector[1] = (size_t)ceil(rhoStep*suppressionVector[1]);
	// 
	for(int count=0;count<(int)linesMax;count++)
	{
		unsigned int maxV = 0;
		size_t maxVi,maxVj;
		// searching for maximum values within the specified range
		if (pre_subs_right<pre_subs_left)
		{
			for (size_t j=0; j<pre_subs_right; j++)
			{
				for (size_t i=0; i<m; i++)
				{
					if (vote[j*m+i]>maxV)
					{
						maxV = vote[j*m+i];
						maxVi = i;
						maxVj = j;
					}
				}
			}
			for (size_t j=pre_subs_left; j<n; j++)
			{
				for (size_t i=0; i<m; i++)
				{
					if (vote[j*m+i]>maxV)
					{
						maxV = vote[j*m+i];
						maxVi = i;
						maxVj = j;
					}
				}
			}
		}
		else
		{
			for (size_t j=pre_subs_left; j<pre_subs_right; j++)
			{
				for (size_t i=0; i<m; i++)
				{
					if (vote[j*m+i]>maxV)
					{
						maxV = vote[j*m+i];
						maxVi = i;
						maxVj = j;
					}
				}
			}
		}
		peaksInfo[count*3] = (double)maxV;
		peaksInfo[count*3+1] = thetaSeq[maxVj];
		peaksInfo[count*3+2] = rhoSeq[maxVi];
		size_t startJ,endJ,startI,endI;
		startJ = 0;
		if (maxVj>suppVector[0])
		{
			startJ = maxVj-suppVector[0];
		}
		endJ = n-1;
		if (maxVj+suppVector[0]<n-1)
		{
			endJ = maxVj+suppVector[0];
		}
		startI = 0;
		if (maxVi>suppVector[1])
		{
			startI = maxVi-suppVector[1];
		}
		endI = 0;
		if (maxVi+suppVector[1]<m-1)
		{
			endI = maxVi+suppVector[1];
		}
		for (size_t j=startJ; j<=endJ; j++)
		{
			for (size_t i=startI; i<=endI; i++)
			{
				vote[j*m+i] = 0;
			}
		}
	}
}
/* Procedure of searching lines */
void lines(unsigned int *img,double *peaksInfo, unsigned int linesMax, unsigned int lineGap, size_t m, size_t n, double *linesInfo)
{
	const int shift = 16;
	for (unsigned int count=0;count<linesMax;count++)
	{
		double theta =  peaksInfo[count*3+1];
		double rho =  peaksInfo[count*3+2];
		linesInfo[count*6+4] = theta;
		linesInfo[count*6+5] = rho;
		double thisSin = sin(theta*pi/180);
		double thisCos = cos(theta*pi/180);
		int rStart = 0;
		int cStart = (int)floor(rho/thisCos);
		// determine whether the cStart overflows the boundary
		if (cStart<0)
		{
			cStart = 0;
			rStart = (int)floor(rho/thisSin);
		}
		else if (cStart>=n)
		{
			cStart = (int)n-1;
			rStart = (int)floor((rho-cStart*thisCos)/thisSin);
		}
		int biasFlag = 0;
		int r0=rStart,c0=cStart,dr0,dc0;
		if (abs(theta)<=45)
		{
			biasFlag = 1;
			dr0 = 1;
			dc0 = (int)floor(fabs(thisSin)*(1 << shift)/thisCos+0.5);
			c0 = (c0 << shift)+(1 << (shift-1));
		}
		else
		{
			dc0 = 1;
			dr0 = (int)floor(thisCos*(1 << shift)/fabs(thisSin)+0.5);
			r0 = (r0 << shift)+(1 << (shift-1));
		}
		if (theta>0)
		{
			dc0 = -dc0;
		}
		// walk along the line using fixed-point arithmetics,
		// ... stop at the image border
		int lastLength = 0, lineStart[2], lineEnd[2];
		for (int gap=0, c=c0, r=r0, startPFlag =1; ; c+=dc0, r+=dr0)
		{
			int r1,c1;
			if (biasFlag)
			{
				c1 = c >> shift;
				r1 = r;
			}
			else
			{
				c1 = c;
				r1 = r >> shift;
			}
			if (c1 < 0||c1 >= n||r1 < 0||r1 >= m)
			{
				// ensure lastLength has been update before  exit
				int thisLength = (int)floor(sqrt(pow((double)lineStart[0]-(double)lineEnd[0],2) + pow((double)lineStart[1]-(double)lineEnd[1],2)));
				if (thisLength > lastLength)
				{
					lastLength = thisLength;
					linesInfo[count*6+0] = (double)lineStart[0];
					linesInfo[count*6+1] = (double)lineStart[1];
					linesInfo[count*6+2] = (double)lineEnd[0];
					linesInfo[count*6+3] = (double)lineEnd[1];
				}
				break;
			}
			if (img[c1*m+r1])
			{
				gap = 0;
				if (startPFlag)
				{
					startPFlag = 0;
					lineStart[0] = r1+1;
					lineStart[1] = c1+1;
					lineEnd[0] = r1+1;
					lineEnd[1] = c1+1;
				}
				else
				{
					lineEnd[0] = r1+1;
					lineEnd[1] = c1+1;
				}
			}
			else if (!startPFlag)
			{
				// fuzzy processing the points beside line
				int leftPFlag = 0, rightPFlag = 0;
				if (c1>0 && img[(c1-1)*m+r1])
					leftPFlag = 1;
				if (c1<n-1 && img[(c1+1)*m+r1])
					rightPFlag = 1;
				if (leftPFlag || rightPFlag)
				{
					gap = 0;
					lineEnd[0] = r1+1;
					lineEnd[1] = c1+1;
					continue;
				}
				// end of line segment
				if (++gap > (int)lineGap)
				{
					startPFlag = 1;
					int thisLength = (int)floor(sqrt(pow((double)lineStart[0]-(double)lineEnd[0],2) + pow((double)lineStart[1]-(double)lineEnd[1],2)));
					if (thisLength > lastLength)
					{
						lastLength = thisLength;
						linesInfo[count*6+0] = (double)lineStart[0];
						linesInfo[count*6+1] = (double)lineStart[1];
						linesInfo[count*6+2] = (double)lineEnd[0];
						linesInfo[count*6+3] = (double)lineEnd[1];
					}
				}
			}
		}
		// enhanced the robustness: none line segment are detected
		if (!lastLength)
		{
			linesInfo[count*6+0] = (double)rStart+1;
			linesInfo[count*6+1] = 1;
			linesInfo[count*6+2] = (double)rStart+1;
			linesInfo[count*6+3] = (double)n;
		}
	}
}
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
// --> Parameters check

    // check for proper number of arguments
    if(nrhs!=8) 
	{
        mexErrMsgIdAndTxt("MyToolbox:fixedsearchLines:nrhs","Eight inputs required.");
    }
    if(nlhs!=1) 
	{
        mexErrMsgIdAndTxt("MyToolbox:fixedsearchLines:nlhs","One output required.");
    }
    // check the type of input argument(s)
    if( !mxIsUint32(prhs[0]) ) 
	{
        mexErrMsgIdAndTxt("MyToolbox:fixedsearchLines:notUint32","Input matrix of image must be type Uint32.");
    }
    if( !mxIsUint32(prhs[1]) ) 
	{
        mexErrMsgIdAndTxt("MyToolbox:fixedsearchLines:notUint32","Input matrix of voting space must be type Uint32.");
    }
    // check the shape of input argument(s)
    if(mxGetNumberOfDimensions(prhs[0])!=2) 
	{
        mexErrMsgIdAndTxt("MyToolbox:fixedsearchLines:not2D","Input of image must be 2-D.");
    }
    if(mxGetNumberOfDimensions(prhs[1])!=2) 
	{
        mexErrMsgIdAndTxt("MyToolbox:fixedsearchLines:not2D","Input of voting space must be 2-D.");
    }
// --> Accept input parameters

    // create a pointer to the real data in the input matrix
    unsigned int *img = mxGetData(prhs[0]);
	unsigned int *vote = mxGetData(prhs[1]);
    // and get dimensions of the input matrix 
	size_t mrowsImg = mxGetM(prhs[0]);
    size_t ncolsImg = mxGetN(prhs[0]);
	size_t mrowsVote = mxGetM(prhs[1]);
    size_t ncolsVote = mxGetN(prhs[1]);
	// other inputs
	double *thetaSeq = mxGetData(prhs[2]);
	double *rhoSeq = mxGetData(prhs[3]);
	unsigned int linesMax = (unsigned int)mxGetScalar(prhs[4]);
	unsigned int lineGap = (unsigned int)mxGetScalar(prhs[5]);
	double angleS = mxGetScalar(prhs[6]);
	double angleE = mxGetScalar(prhs[7]);
    // check the range of input argument(s)
	if (angleS<-90 || angleS>90 || angleE<-90 || angleE>90)
	{
        mexErrMsgIdAndTxt("MyToolbox:fixedsearchLines:OutOfDefined","Input of specified range must be between [-90,90].");
	}
// --> Parameters can be added to function interface

	// default parameters
	
	// check them
	
// --> Prepare data
	
	// prepare information of step of theta & rho
	double thetaStep = thetaSeq[1]-thetaSeq[0];
	double rhoStep = rhoSeq[1]-rhoSeq[0];
	size_t thetaQ = 45;
	// creat the matrix to store information of peaks
	double *peaksInfo = mxMalloc(3*linesMax*sizeof(double));
	// creat the matrix to store information of lines
    plhs[0] = mxCreateDoubleMatrix (6,linesMax,mxREAL);
    double *linesInfo = mxGetData(plhs[0]);
	// compute position of angleS and angleE
	int count =0;
	for (int base = -45;angleS>=base;count++)
	{
		base+=45;
	}
	size_t angleSQ = (size_t)count*thetaQ+(size_t)(floor((angleS+(double)(2-count)*45)/thetaStep));
	count = 0;
	for (int base = -45;angleE>=base;count++)
	{
		base+=45;
	}
	size_t angleEQ = (size_t)count*thetaQ+(size_t)(ceil((angleE+(double)(2-count)*45)/thetaStep));
	if (angleEQ>4*thetaQ)
	{
		angleEQ-=4*thetaQ;
	}
// --> Main computation

    // call the procedure of searching peaks
    peaks(vote,mrowsVote,ncolsVote,linesMax,thetaSeq,rhoSeq,thetaStep,rhoStep,peaksInfo,angleSQ,angleEQ);
	// call the procedure of searching lines
	lines(img,peaksInfo,linesMax,lineGap,mrowsImg,ncolsImg,linesInfo);
// --> Rest

	// release memory
	mxFree(peaksInfo);
}