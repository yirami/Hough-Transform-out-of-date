/* ===============================================================
 * Copyright (C) 2017, Yirami  
 * Released under BSD License
 * All rights reserved.
 * 
 * Filename: votePPHT.c (compiled into mex)
 * F-Description: Progressive Probabilistic Hough Transform(PPHT)
 *					with a limit of angle
 *
 *		inputs 	a MxN matrix "img"
 *				a 	  scalar "linesMax"
 *				a 	  scalar "lineLength"
 *				a 	  scalar "lineGap"
 *			 	a 	  scalar "angleS"
 *			 	a	  scalar "angleE"
 * 			and 
 *		outputs a 4xR matrix "lines"
 *
 * 		The calling syntax is:
 *			[linesEnd] = votePPHTa(img,linesMax,lineLength,lineGap,angleS,angleE)
 * 
 * Author：Yirami
 * A-Description:
 * 		mailto:	  1446675043@qq.com
 * 		website:  https://yirami.github.io./
 * 
 * Date: 24-05-2017
 * Version: 1.0
 *  
 * New  Release:
 * N-Description:
 *
 * =============================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include <stdlib.h>
#include <time.h>

#define pi 3.14159265358979

/* Traversal existing endpoints */
unsigned int traversalLoop (unsigned int confusion, unsigned int linesCountCurr, double *linesPolar, double thetaCurr, double rhoCurr, unsigned int *hit)
{
	unsigned int hitFlag = 0;
	for (unsigned int count = 0; count <linesCountCurr;count++)
	{
		/*double diffTR = 0.5*pow(linesPolar[2*count]-thetaCurr,2)+0.2*pow(linesPolar[2*count+1]-rhoCurr,2);
		if (diffTR <= (double)confusion)
		{
			hitFlag = 1;
			*hit = count;
			break;
		}*/
		if ( abs(linesPolar[2*count]-thetaCurr) < 0.5*(double)confusion && abs(linesPolar[2*count+1]-rhoCurr) <= (double)confusion )
		{
			hitFlag = 1;
			*hit = count;
			break;
		}
	}
	return hitFlag;
}

/* The computational routine */
void votePPHTa(	unsigned int *img, unsigned int *vote, size_t m, size_t n, \
				double *sinMap, double *cosMap, double thetaStep, \
				unsigned int thetaQ, double rhoStep, unsigned int rhoQ, \
				unsigned int *outputInfo[2], unsigned int lineLimitInfo[6], \
				unsigned int pre_subs_left, unsigned int pre_subs_right	)
{
// --> Unpackage the input parameters

	// lineLimitInfo[6]
	unsigned int threshold = lineLimitInfo[0];
	int lineLength = (int)lineLimitInfo[1];
	unsigned int lineGap = lineLimitInfo[2];
	unsigned int linesMax = lineLimitInfo[3];
	unsigned int isMerge = lineLimitInfo[4];
	unsigned int confusion = lineLimitInfo[5];
	if (isMerge)
		double *linesPolar = mxMalloc(2*linesMax*sizeof(double));
	// outputInfo[2]
	unsigned int *linesCount = outputInfo[0];
	unsigned int *linesEndSrc = outputInfo[1];
// --> Stage 1. collect non-zero image points

	int edgePointsCount = 0, defaultSize = 1000;
	unsigned int *edgePoints = mxMalloc(2*defaultSize*sizeof(unsigned int));
	// init "mask" : 
	//				1 means doesn't delong to any line 
	// 			and 
	// 				2 means has been voted
	int *mask = mxCalloc(m*n,sizeof(int));
	// "j" means cols-1
    for (unsigned int j=0; j<n; j++) 
	{
		// and "i" means rows-1
		for (unsigned int i=0; i<m; i++)
		{
			if (img[j*m+i])
			{
				mask[j*m+i] = 1;
				if (edgePointsCount%defaultSize==0&&edgePointsCount!=0)
				{
					edgePoints = mxRealloc(edgePoints, 2*(edgePointsCount/defaultSize+1)*defaultSize*sizeof(unsigned int));
				}
				edgePoints[2*edgePointsCount] = i;
				edgePoints[2*edgePointsCount+1] = j;
				edgePointsCount++;
			}
		}
    }
// --> Stage 2. process all the points in random order

	for (int count = edgePointsCount; count>0; count--)
	{
		// choose random point out of the remaining ones
		srand((unsigned)time(NULL));
		int index = rand()%count;
		unsigned int voteMax = threshold-1, numMax;
		double thetaCurr, rhoCurr;
		unsigned int point[2] = {edgePoints[2*index],edgePoints[2*index+1]};
		int lineEnd[4] = {0};
		// "remove" it by overriding it with the last element
		edgePoints[2*index] = edgePoints[2*(count-1)];
		edgePoints[2*index+1] = edgePoints[2*(count-1)+1];
		// check if it has been excluded already (i.e. belongs to some other line)
		if (!mask[m*point[1]+point[0]])
			continue;
		// update accumulator, find the most probable line
		// ... with a limit of angle
		if (pre_subs_right	<pre_subs_left)
		{
			unsigned int val_curr;
			for (unsigned int idx = 0; idx < pre_subs_right; idx++)
			{
				double rho_c = floor((point[1]*cosMap[idx]+point[0]*sinMap[idx])/rhoStep+0.5);
				val_curr = ++vote[idx*(2*rhoQ+1)+(unsigned int)rho_c+rhoQ];
				if (voteMax < val_curr)
				{
					voteMax = val_curr;
					numMax = idx;
					if (isMerge)
					{
						thetaCurr = -90+45*(idx/thetaQ)+thetaStep*(idx%thetaQ);
						rhoCurr = rho_c;
					}
				}
			}
			for (unsigned int idx = pre_subs_left; idx < 4*thetaQ; idx++)
			{
				double rho_c = floor((point[1]*cosMap[idx]+point[0]*sinMap[idx])/rhoStep+0.5);
				val_curr = ++vote[idx*(2*rhoQ+1)+(unsigned int)rho_c+rhoQ];
				if (voteMax < val_curr)
				{
					voteMax = val_curr;
					numMax = idx;
					if (isMerge)
					{
						thetaCurr = -90+45*(idx/thetaQ)+thetaStep*(idx%thetaQ);
						rhoCurr = rho_c;
					}
				}
			}
		}
		else
		{
			for (unsigned int idx = pre_subs_left; idx < pre_subs_right; idx++)
			{
				double rho_c = floor((point[1]*cosMap[idx]+point[0]*sinMap[idx])/rhoStep+0.5);
				unsigned int val_curr = ++vote[idx*(2*rhoQ+1)+(unsigned int)rho_c+rhoQ];
				if (voteMax < val_curr)
				{
					voteMax = val_curr;
					numMax = idx;
					if (isMerge)
					{
						thetaCurr = -90+45*(idx/thetaQ)+thetaStep*(idx%thetaQ);
						rhoCurr = rho_c;
					}
				}
			}
		}
		// record here has been vote
		mask[m*point[1]+point[0]]++;		
		// if it is too "weak" candidate, continue with another point
		if (voteMax < threshold)
			continue;
		// from the current point walk in each direction 
		// ... along the found line and extract the line segment
        int biasFlag, r0 = point[0], c0 = point[1], dr0, dc0;
		const int shift = 16;
		double a = sinMap[numMax], b = cosMap[numMax];
		// if "dc > dr"
		if (fabs(a) > fabs(b))
		{
			biasFlag = 1;
			// ensure dr0 > 0
			dc0 = a < 0 ? 1 : -1;
			dr0 = (int)floor(b*(1 << shift)/fabs(a)+0.5);
			r0 = (r0 << shift)+(1 << (shift-1));
		}
		// if "dr >= dc"
		else
		{
			biasFlag = 0;
			// ensure dc0 > 0
			dr0 = a < 0 ? 1 : -1;
			dc0 = (int)floor(fabs(a)*(1 << shift)/b+0.5);
			c0 = (c0 << shift)+(1 << (shift-1));
		}
		// search the end of this line
		for (int k = 0; k < 2; k++)
		{
			unsigned int gap = 0;
			int c =c0, r = r0, dc = dc0, dr = dr0;
			if (k > 0)
				dc = -dc, dr = -dr;
			// walk along the line using fixed-point arithmetics,
			// ... stop at the image border or in case of too big gap
			for (;; c += dc, r += dr)
			{
				int r1,c1;
				if (biasFlag)
				{
					c1 = c;
					r1 = r >> shift;
				}
				else
				{
					c1 = c >> shift;
					r1 = r;
				}
				if (c1 < 0||c1 >= n||r1 < 0||r1 >= m)
					break;
				// for each point meet the requirement :
				// 		clear the gap;
				//		record the coordinate;
				if (mask[c1*m+r1])
				{
					gap = 0;
					lineEnd[2*k] = r1;
					lineEnd[2*k+1] = c1;
				}
				else if (++gap > lineGap)
					break;
			}
		}
		// determine whether the line meets the requirements
		int goodLine = 0;
		if (abs(lineEnd[0]-lineEnd[2]) >= lineLength||abs(lineEnd[1]-lineEnd[3]) >= lineLength)
			goodLine = 1;
		// traverse all points on the line and perform appropriate actions
		for (int k = 0; k < 2; k++)
		{
			int c =c0, r = r0, dc = dc0, dr = dr0;
			if (k > 0)
				dc = -dc, dr = -dr;
			// walk along the line using fixed-point arithmetics,
			// ... stop at the endpoint
			for (;; c += dc, r += dr)
			{
				int r1,c1;
				if (biasFlag)
				{
					c1 = c;
					r1 = r >> shift;
				}
				else
				{
					c1 = c >> shift;
					r1 = r;
				}
				if (c1 < 0||c1 >= n||r1 < 0||r1 >= m)
					break;
				// for each point meet the requirement :
				// 		unvote;
				//		clear the "mask";
				if (mask[c1*m+r1])
				{
					if (goodLine && mask[c1*m+r1]==2)
					{

						if (pre_subs_right<pre_subs_left)
						{
							for (unsigned int idx = 0; idx < pre_subs_right; idx++)
							{
								double rho_c = floor((c1*cosMap[idx]+r1*sinMap[idx])/rhoStep+0.5);
								vote[idx*(2*rhoQ+1)+(unsigned int)rho_c+rhoQ]--;
							}
							for (unsigned int idx = pre_subs_left; idx < 4*thetaQ; idx++)
							{
								double rho_c = floor((c1*cosMap[idx]+r1*sinMap[idx])/rhoStep+0.5);
								vote[idx*(2*rhoQ+1)+(unsigned int)rho_c+rhoQ]--;
							}
						}
						else
						{
							for (unsigned int idx = pre_subs_left; idx < pre_subs_right; idx++)
							{
								double rho_c = floor((c1*cosMap[idx]+r1*sinMap[idx])/rhoStep+0.5);
								vote[idx*(2*rhoQ+1)+(unsigned int)rho_c+rhoQ]--;
							}
						}
					}
					mask[c1*m+r1] = 0;
				}
				if (r1 == lineEnd[2*k] && c1 == lineEnd[2*k+1])
					break;
			}
		}
		// record the endpoint
		if (goodLine)
		{
			if (isMerge)
			{
				unsigned int hitNum = 0;
				if (traversalLoop (confusion,*linesCount,linesPolar,thetaCurr,rhoCurr,&hitNum))
				{
					// search the most side endpoint
					unsigned int leftP[2] = {n,n}, rightP[2] = {0};
					
					if (linesEndSrc[4*hitNum+1] < leftP[1])
					{
						leftP[0] = linesEndSrc[4*hitNum];
						leftP[1] = linesEndSrc[4*hitNum+1];
					}
					if (linesEndSrc[4*hitNum+3] < leftP[1])
					{
						leftP[0] = linesEndSrc[4*hitNum+2];
						leftP[1] = linesEndSrc[4*hitNum+3];
					}
					if ((unsigned int)lineEnd[1] < leftP[1])
					{
						leftP[0] = (unsigned int)lineEnd[0];
						leftP[1] = (unsigned int)lineEnd[1];
					}
					if ((unsigned int)lineEnd[3] < leftP[1])
					{
						leftP[0] = (unsigned int)lineEnd[2];
						leftP[1] = (unsigned int)lineEnd[3];
					}
					if (linesEndSrc[4*hitNum+1] > rightP[1])
					{
						rightP[0] = linesEndSrc[4*hitNum];
						rightP[1] = linesEndSrc[4*hitNum+1];
					}
					if (linesEndSrc[4*hitNum+3] > rightP[1])
					{
						rightP[0] = linesEndSrc[4*hitNum+2];
						rightP[1] = linesEndSrc[4*hitNum+3];
					}
					if ((unsigned int)lineEnd[1] > rightP[1])
					{
						rightP[0] = (unsigned int)lineEnd[0];
						rightP[1] = (unsigned int)lineEnd[1];
					}
					if ((unsigned int)lineEnd[3] > rightP[1])
					{
						rightP[0] = (unsigned int)lineEnd[2];
						rightP[1] = (unsigned int)lineEnd[3];
					}					
					// replace the data
					linesEndSrc[4*hitNum] = leftP[0];
					linesEndSrc[4*hitNum+1] = leftP[1];
					linesEndSrc[4*hitNum+2] = rightP[0];
					linesEndSrc[4*hitNum+3] = rightP[1];
					double angle = atan(((double)rightP[1]-(double)leftP[1])/((double)leftP[0]-(double)rightP[0]));
					linesPolar[2*hitNum] = angle;
					linesPolar[2*hitNum+1] = floor((leftP[1]*cos(angle)+leftP[0]*sin(angle))/rhoStep+0.5);
				}
				else
				{
					linesEndSrc[4*(*linesCount)] = (unsigned int)lineEnd[0];
					linesEndSrc[4*(*linesCount)+1] = (unsigned int)lineEnd[1];
					linesEndSrc[4*(*linesCount)+2] = (unsigned int)lineEnd[2];
					linesEndSrc[4*(*linesCount)+3] = (unsigned int)lineEnd[3];
					linesPolar[2*(*linesCount)] = thetaCurr;
					linesPolar[2*(*linesCount)+1] = rhoCurr;
					(*linesCount)++;
					if ( *linesCount >= linesMax )
					{
						mxFree(edgePoints);
						mxFree(mask);
						mxFree(linesPolar);
						return;
					}
				}
			}
			else
			{
				linesEndSrc[4*(*linesCount)] = (unsigned int)lineEnd[0];
				linesEndSrc[4*(*linesCount)+1] = (unsigned int)lineEnd[1];
				linesEndSrc[4*(*linesCount)+2] = (unsigned int)lineEnd[2];
				linesEndSrc[4*(*linesCount)+3] = (unsigned int)lineEnd[3];
				(*linesCount)++;
				if ( *linesCount >= linesMax )
				{
					mxFree(edgePoints);
					mxFree(mask);
					return;
				}
			}
		}
	}
	// release dynamically allocated memory 
	mxFree(edgePoints);
	mxFree(mask);
	if (isMerge)
		mxFree(linesPolar);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
// --> Parameters check

    // check for proper number of arguments
    if(nrhs!=6) 
	{
        mexErrMsgIdAndTxt("MyToolbox:votePPHT:nrhs","Six inputs required.");
    }
    if(nlhs!=1) 
	{
        mexErrMsgIdAndTxt("MyToolbox:votePPHT:nlhs","One output required.");
    }
    // check the type of input argument(s)
    if( !mxIsUint32(prhs[0]) ) 
	{
        mexErrMsgIdAndTxt("MyToolbox:votePPHT:notUint32","Input matrix must be type Uint32.");
    }
    // check the shape of input argument(s)
    if(mxGetNumberOfDimensions(prhs[0])!=2) 
	{
        mexErrMsgIdAndTxt("MyToolbox:votePPHT:not2D","Input must be 2-D.");
    }
// --> Accept input parameters

    // create a pointer to the real data in the input matrix
    unsigned int *img = mxGetData(prhs[0]);
    // and get dimensions of the input matrix 
	size_t mrows = mxGetM(prhs[0]);
    size_t ncols = mxGetN(prhs[0]);
	// add to
	unsigned int linesMax = (unsigned int)mxGetScalar(prhs[1]);
	unsigned int lineLength = (unsigned int)mxGetScalar(prhs[2]);
	unsigned int lineGap = (unsigned int)mxGetScalar(prhs[3]);
    // get the upper and lower of angle
	double angleS = mxGetScalar(prhs[4]);
	double angleE = mxGetScalar(prhs[5]);
// --> Parameters can be added to function interface

	// default parameters
	double thetaResolution = 1;
	double rhoResolution = 1;
	unsigned int threshold = 20;
	//unsigned int lineLength = 200;
	//unsigned int lineGap = 5;
	//unsigned int linesMax = 100;
	unsigned int isMerge = 0;
	unsigned int confusion = 4;
	// check them
	double rhoMax = floor(sqrt(pow((double)mrows-1,2)+pow((double)ncols-1,2))+0.5);
    if(thetaResolution>5||rhoResolution>rhoMax/10) 
	{
        mexErrMsgIdAndTxt("MyToolbox:votePPHT:meaningless","Resolution of Theta or Rho is meaningless.");
    }
// --> 	Prepare data
	
	// compute quantity and step of "theta" & "rho"
	unsigned int thetaQ = (unsigned int)ceil(45/thetaResolution);
	unsigned int rhoQ = (unsigned int)ceil(rhoMax/rhoResolution);
	double thetaStep = (double)45/thetaQ;
	double rhoStep = (double)rhoMax/rhoQ;
	// create "H" & "theta" & "rho" and pass the pointer
    mxArray *HPtr = mxCreateNumericMatrix(2*rhoQ+1,4*thetaQ,mxUINT32_CLASS,mxREAL);
    mxArray *thetaPtr = mxCreateDoubleMatrix (1,4*thetaQ,mxREAL);
    mxArray *rhoPtr = mxCreateDoubleMatrix (1,2*rhoQ+1,mxREAL);
	unsigned int *H = mxGetData(HPtr);
	double *theta = mxGetPr(thetaPtr);
	double *rho = mxGetPr(rhoPtr);
	// dynamic memory allocation of "sinMap" and "cosMap"
    double *sinMap = mxMalloc(4*thetaQ*sizeof(double));
    double *cosMap = mxMalloc(4*thetaQ*sizeof(double));
	// compute "theta" & "sinMap" & "cosMap"
	for (unsigned int count = 0; count<thetaQ; count++)
	{
		unsigned int start[4] = {0, thetaQ, 2*thetaQ, 3*thetaQ};
		// compute theta
		theta[start[0]+count] = -90 + thetaStep*count;
		theta[start[1]+count] = -45 + thetaStep*count;
		theta[start[2]+count] =  0  + thetaStep*count;
		theta[start[3]+count] =  45 + thetaStep*count;
		// compute sinMap
		sinMap[start[0]+count] = sin(theta[start[0]+count]*pi/180);
		sinMap[start[1]+count] = sin(theta[start[1]+count]*pi/180);
		sinMap[start[2]+count] = sin(theta[start[2]+count]*pi/180);
		sinMap[start[3]+count] = sin(theta[start[3]+count]*pi/180);
		// compute cosMap
		cosMap[start[0]+count] = cos(theta[start[0]+count]*pi/180);
		cosMap[start[1]+count] = cos(theta[start[1]+count]*pi/180);
		cosMap[start[2]+count] = cos(theta[start[2]+count]*pi/180);
		cosMap[start[3]+count] = cos(theta[start[3]+count]*pi/180);
	}
    // compute  "rho" 
	for (unsigned int count = 1; count<=rhoQ; count++)
	{
		rho[rhoQ+count] = rhoStep*(double)count;
		rho[rhoQ-count] = -rhoStep*(double)count;
	}
	// dynamic allocation of "linesEnd" & ""
	unsigned int *linesEndSrc = mxMalloc(4*linesMax*sizeof(unsigned int));
	// compute position of angleS and angleE
	unsigned int count =0;
	for (int base = -45;angleS>=base;count++)
	{
		base+=45;
	}
	unsigned int angleSQ = count*thetaQ+(unsigned int)(floor((angleS+(double)(2-count)*45)/thetaStep));
	count = 0;
	for (int base = -45;angleE>=base;count++)
	{
		base+=45;
	}
	unsigned int angleEQ = count*thetaQ+(unsigned int)(ceil((angleE+(double)(2-count)*45)/thetaStep));
	if (angleEQ>4*thetaQ)
	{
		angleEQ-=4*thetaQ;
	}
// --> Main computation

	// package the input parameters
	unsigned int lineLimitInfo[6] = 
	{
		threshold,
		lineLength,
		lineGap,
		linesMax,
		isMerge,
		confusion
	};
    unsigned int linesNumber = 0;
	unsigned int *outputInfo[2] =
	{
		&linesNumber,
		linesEndSrc
	};
    // call the computational routine
    votePPHTa(img,H,mrows,ncols,sinMap,cosMap,thetaStep,thetaQ,rhoStep,rhoQ,outputInfo,lineLimitInfo,angleSQ,angleEQ);
// --> Rest

	// create the output matrix "linesEnd"
	plhs[0] = mxCreateNumericMatrix(4,linesNumber,mxUINT32_CLASS,mxREAL);
    unsigned int *linesEnd = mxGetData(plhs[0]);
	// and copy the source by pointer
	memcpy(linesEnd, linesEndSrc, 4*linesNumber*sizeof(unsigned int));
	// release memory
	mxFree(sinMap);
	mxFree(cosMap);
	mxFree(linesEndSrc);
	mxDestroyArray(HPtr);
	mxDestroyArray(thetaPtr);
	mxDestroyArray(rhoPtr);
}
