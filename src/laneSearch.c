/* ===============================================================
 * Copyright (C) 2017, Yirami  
 * Released under BSD License
 * All rights reserved.
 * 
 * Filename: laneSearch.c (compiled into mex)
 * F-Description: search lane in global scope
 *
 *		inputs 	a MxN matrix "img"
 *			 	a PxQ matrix "vote"
 *			 	a 	  vector "thetaSeq"
 *			 	a 	  vector "rhoSeq"
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
 *			[linesInfo] = laneSearch(img,vote,thetaSeq,rhoSeq)
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
#define batchQ 9
#define timesEachBatch 4

/* Procedure of clean */
void cleanVote(unsigned int *vote, size_t rCenter, size_t cCenter, size_t *areaSpec, size_t rWhole, size_t cWhole)
{
	// determine the legal start and end of rows and column
	size_t startJ,endJ,startI,endI;
	startJ = 0;
	if (cCenter>areaSpec[0])
	{
		startJ = cCenter-areaSpec[0];
	}
	endJ = cWhole-1;
	if (cCenter+areaSpec[0]<cWhole-1)
	{
		endJ = cCenter+areaSpec[0];
	}
	startI = 0;
	if (rCenter>areaSpec[1])
	{
		startI = rCenter-areaSpec[1];
	}
	endI = rWhole-1;
	if (rCenter+areaSpec[1]<rWhole-1)
	{
		endI = rCenter+areaSpec[1];
	}
	// clean the specific area
	for (size_t j=startJ; j<=endJ; j++)
	{
		for (size_t i=startI; i<=endI; i++)
		{
			vote[j*rWhole+i] = 0;
		}
	}
}

/* Procedure of searching peaks */
void peaks(unsigned int *vote, size_t m, size_t n, double *thetaSeq, double*rhoSeq, double thetaStep, double rhoStep, double *peaksInfo, size_t *actualPeaksNumPtr)
{
	// built-in parameters
	double suppressionVector[2] = {2,4};
	int threshold = 50;
	unsigned int batchPriorityVector[batchQ] = {5,4,6,3,7,2,8,1,9};
		// Note: 1 -> -90~-70; 5 -> -10~10 ...
	unsigned int breakpointBatch = 4;
		// Note: in here can make a decision to quit
	// prepare data
	size_t suppVector[2];
	suppVector[0] = (size_t)ceil(thetaStep*suppressionVector[0]);
	suppVector[1] = (size_t)ceil(rhoStep*suppressionVector[1]);
	unsigned int batchSize = (unsigned int)floor((double)n/(double)batchQ);
	unsigned int batchBias = (unsigned int)floor((n-batchQ*batchSize)/2);
	// batch loop
	for (int batch=0;batch<batchQ;batch++)
	{
		unsigned int thisSeq = batchPriorityVector[batch]-1;
		unsigned int base = batchBias+thisSeq*batchSize;
		double batchPeaks[3*timesEachBatch];
		size_t batchPeaksNum = 0;
		// repeat batch searching
		for (int count=0;count<timesEachBatch;count++)
		{
			unsigned int maxV = threshold;
			size_t maxVi,maxVj;
			unsigned int hitFlag =0;
			for (size_t j=base; j<base+batchSize; j++)
			{
				for (size_t i=0; i<m; i++)
				{
					if (vote[j*m+i]>maxV)
					{
						hitFlag = 1;
						maxV = vote[j*m+i];
						maxVi = i;
						maxVj = j;
					}
				}
			}
			// if can't find at this time then quit this batch
			if (!hitFlag)
			{
				break;
			}
			// else store the infomation of maximum point
			batchPeaksNum++;
			batchPeaks[3*count] = (double)maxV;
			batchPeaks[3*count+1] = thetaSeq[maxVj];
			batchPeaks[3*count+2] = rhoSeq[maxVi];
			// clear the vote space near this maximum point
			cleanVote(vote,maxVi,maxVj,&suppVector,m,n);
		}
		// output current results
		for (int count=0;count<batchPeaksNum;count++)
		{
			peaksInfo[3*actualPeaksNumPtr[0]] = batchPeaks[3*count];
			peaksInfo[3*actualPeaksNumPtr[0]+1] = batchPeaks[3*count+1];
			peaksInfo[3*actualPeaksNumPtr[0]+2] = batchPeaks[3*count+2];
			actualPeaksNumPtr[0]++;
		}
		// deal with breakpoint
		if (batch==breakpointBatch)
		{
			break;
		}
	}
}
/* Procedure of searching lines */
void lines(unsigned int *img,double *peaksInfo, size_t actualPeaksNum, unsigned int lineGap, size_t m, size_t n, double *linesInfo, size_t *actualLinesNumPtr)
{
	const int shift = 16;
	for (unsigned int count=0;count<actualPeaksNum;count++)
	{
		// load the theta anf rho
		double theta =  peaksInfo[count*3+1];
		double rho =  peaksInfo[count*3+2];
		// output the theta anf rho
		linesInfo[actualLinesNumPtr[0]*6+4] = theta;
		linesInfo[actualLinesNumPtr[0]*6+5] = rho;
		// prepare data
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
					linesInfo[actualLinesNumPtr[0]*6+0] = (double)lineStart[0]+1;
					linesInfo[actualLinesNumPtr[0]*6+1] = (double)lineStart[1]+1;
					linesInfo[actualLinesNumPtr[0]*6+2] = (double)lineEnd[0]+1;
					linesInfo[actualLinesNumPtr[0]*6+3] = (double)lineEnd[1]+1;
				}
				break;
			}
			if (img[c1*m+r1])
			{
				gap = 0;
				if (startPFlag)
				{
					startPFlag = 0;
					lineStart[0] = r1;
					lineStart[1] = c1;
					lineEnd[0] = r1;
					lineEnd[1] = c1;
				}
				else
				{
					lineEnd[0] = r1;
					lineEnd[1] = c1;
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
					lineEnd[0] = r1;
					lineEnd[1] = c1;
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
						linesInfo[actualLinesNumPtr[0]*6+0] = (double)lineStart[0]+1;
						linesInfo[actualLinesNumPtr[0]*6+1] = (double)lineStart[1]+1;
						linesInfo[actualLinesNumPtr[0]*6+2] = (double)lineEnd[0]+1;
						linesInfo[actualLinesNumPtr[0]*6+3] = (double)lineEnd[1]+1;
					}
				}
			}
		}
		// enhanced the robustness: none line segment are detected
		if (lastLength)
		{
			actualLinesNumPtr[0]++;
		}
	}
}
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
// --> Parameters check

    // check for proper number of arguments
    if(nrhs!=4) 
	{
        mexErrMsgIdAndTxt("MyToolbox:laneSearch:nrhs","Four inputs required.");
    }
    if(nlhs!=1) 
	{
        mexErrMsgIdAndTxt("MyToolbox:laneSearch:nlhs","One output required.");
    }
    // check the type of input argument(s)
    if( !mxIsUint32(prhs[0]) ) 
	{
        mexErrMsgIdAndTxt("MyToolbox:laneSearch:notUint32","Input matrix of image must be type Uint32.");
    }
    if( !mxIsUint32(prhs[1]) ) 
	{
        mexErrMsgIdAndTxt("MyToolbox:laneSearch:notUint32","Input matrix of voting space must be type Uint32.");
    }
    // check the shape of input argument(s)
    if(mxGetNumberOfDimensions(prhs[0])!=2) 
	{
        mexErrMsgIdAndTxt("MyToolbox:laneSearch:not2D","Input of image must be 2-D.");
    }
    if(mxGetNumberOfDimensions(prhs[1])!=2) 
	{
        mexErrMsgIdAndTxt("MyToolbox:laneSearch:not2D","Input of voting space must be 2-D.");
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
// --> Parameters can be added to function interface

	// default parameters
	size_t actualPeaksNum = 0;
	size_t actualLinesNum = 0;
	unsigned int lineGap = 5;
	// check them
	
// --> Prepare data
	
	// prepare information of step of theta & rho
	double thetaStep = thetaSeq[1]-thetaSeq[0];
	double rhoStep = rhoSeq[1]-rhoSeq[0];
	// creat the matrix to store information of peaks
	double *peaksInfo = mxMalloc(3*timesEachBatch*batchQ*sizeof(double));
	
// --> Main computation

    // call the procedure of searching peaks
    peaks(vote,mrowsVote,ncolsVote,thetaSeq,rhoSeq,thetaStep,rhoStep,peaksInfo,&actualPeaksNum);
	// adjust matrix to actual size
	peaksInfo = mxRealloc(peaksInfo,3*actualPeaksNum*sizeof(double));
	// creat the matrix to store information of lines
	double *linesInfo = mxMalloc(6*actualPeaksNum*sizeof(double));
	// call the procedure of searching lines
	lines(img,peaksInfo,actualPeaksNum,lineGap,mrowsImg,ncolsImg,linesInfo,&actualLinesNum);
// --> Rest

	// creat the matrix to store information of lines and make a copy
    plhs[0] = mxCreateDoubleMatrix (6,actualLinesNum,mxREAL);
    double *linesINFO = mxGetData(plhs[0]);
	linesINFO = memcpy(linesINFO,linesInfo,6*actualLinesNum*sizeof(double));

	// release memory
	mxFree(peaksInfo);
	mxFree(linesInfo);
}