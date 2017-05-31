/* ===============================================================
 * Copyright (C) 2017, Yirami  
 * Released under BSD License
 * All rights reserved.
 * 
 * Filename: voteDCHT.c (compiled into mex)
 * F-Description: vote for the algorithm presented in this paper
 *
 *		input 	a MxN matrix "img"
 * 			and 
 *		outputs a PxQ matrix "H"
 *				a 1xR matrix "theta"
 *				a 1xR matrix "rho"
 *
 * 		The calling syntax is:
 *			[H,theta,rho] = voteDCHT(img)
 * 
 * Author£ºYirami
 * A-Description:
 * 		mailto:	  1446675043@qq.com
 * 		website:  https://yirami.github.io./
 * 
 * Date: 29-03-2017
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
void voteDCHT(unsigned int *img, unsigned int *vote, size_t m, size_t n, double rho_step, double *sinMap, double *cosMap, unsigned int thetaQ, unsigned int rhoQ)
{
// --> Vote based on angle range
	// "j" means cols-1
	for (unsigned int j=1; j<(n-1); j++) 
	{
		// and "i" means rows-1
		for (unsigned int i=1;i<(m-1);i++)
		{
			if (img[j*m+i])
			{
				// index of start angle
				unsigned int subs_left;
				// index of end angle
				unsigned int subs_right;
				// flag of skip the process of part of default
				int skipFlag = 1;
				// index of reference point
				unsigned int subs_ref = (j-1)*m+i-1;
				// convolution result
				unsigned int sum = 1000*(img[subs_ref]+img[subs_ref+2*m+2])+ \
									100*(img[subs_ref+m]+img[subs_ref+m+2])+ \
									10*(img[subs_ref+2]+img[subs_ref+2*m])+  \
									img[subs_ref+1]+img[subs_ref+2*m+1];
				// for each point not zero judgment angle range based on convolution
				switch (sum)
				{
					// Note: When we need to reduce the code space, we can compress
					// 	... the range of the "case" to reduce the size of the table !
					//  ... DETAILS may get here : 
					//  ... 	http://www.cnblogs.com/idorax/p/6275259.html
					//
					//  The configuration below can be adjusted as appropriate
					//
					
					// -90~-45
					case 1002:case 2001:case 1001:case 2002:
						subs_left = 0;
						subs_right = thetaQ;
						break;
					// -45~0
					case 2100:case 1200:case 1100:case 2200:
						subs_left = thetaQ;
						subs_right = 2*thetaQ;
						break;
					// 0~45
					case 210:case 120:case 110:case 220:
						subs_left = 2*thetaQ;
						subs_right = 3*thetaQ;
						break;
					// 45~90
					case 21:case 12:case 11:case 22:
						subs_left = 3*thetaQ;
						subs_right = 4*thetaQ;
						break;
					// -90~0
					case 1101:case 2000:case 2101:case 1201:case 1102:case 2111:
						subs_left = 0;
						subs_right = 2*thetaQ;
						break;
					// -45~45
					case 1110:case 200:case 2110:case 1210:case 1120:case 1211:
						subs_left = thetaQ;
						subs_right = 3*thetaQ;
						break;
					// 0~90
					case 111:case 20:case 211:case 121:case 112:case 1121:
						subs_left = 2*thetaQ;
						subs_right = 4*thetaQ;
						break;
					// 45~-45
					case 1011:case 2:case 2011:case 1021:case 1012:case 1112:
						subs_left = 3*thetaQ;
						subs_right = thetaQ;
						break;
					// else
					default:
						if (skipFlag)
						{
							continue;
						}
						else
						{
							subs_left = 0;
							subs_right = 4*thetaQ;
							break;
						}
				}
				// vote based on angle range
				if (subs_right<subs_left)
				{
					for (unsigned int idx = 0; idx < subs_right; idx++)
					{
						double rho_c = floor((j*cosMap[idx]+i*sinMap[idx])/rho_step+0.5);
						vote[idx*(2*rhoQ+1)+(unsigned int)rho_c+rhoQ]++;
					}
					for (unsigned int idx = subs_left; idx < 4*thetaQ; idx++)
					{
						double rho_c = floor((j*cosMap[idx]+i*sinMap[idx])/rho_step+0.5);
						vote[idx*(2*rhoQ+1)+(unsigned int)rho_c+rhoQ]++;
					}
				}
				else
				{
					for (unsigned int idx = subs_left; idx < subs_right; idx++)
					{
						double rho_c = floor((j*cosMap[idx]+i*sinMap[idx])/rho_step+0.5);
						vote[idx*(2*rhoQ+1)+(unsigned int)rho_c+rhoQ]++;
					}
				}
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
        mexErrMsgIdAndTxt("MyToolbox:voteDCHT:nrhs","One input required.");
    }
    if(nlhs!=3) 
	{
        mexErrMsgIdAndTxt("MyToolbox:voteDCHT:nlhs","Three outputs required.");
    }
    // check the type of input argument(s)
    if( !mxIsUint32(prhs[0]) ) 
	{
        mexErrMsgIdAndTxt("MyToolbox:voteDCHT:notUint32","Input matrix must be type Uint32.");
    }
    // check the shape of input argument(s)
    if(mxGetNumberOfDimensions(prhs[0])!=2) 
	{
        mexErrMsgIdAndTxt("MyToolbox:voteDCHT:not2D","Input must be 2-D.");
    }
// --> Accept input parameters

    // create a pointer to the real data in the input matrix
    unsigned int *img = mxGetData(prhs[0]);
    // and get dimensions of the input matrix 
	size_t mrows = mxGetM(prhs[0]);
    size_t ncols = mxGetN(prhs[0]);
    // check the size of matrix are more than 3x3
    if(mrows<3||ncols<3) 
	{
        mexErrMsgIdAndTxt("MyToolbox:voteDCHT:toosmall","Input matrix must be larger than 3x3.");
    }
// --> Parameters can be added to function interface

	// default parameters
	double thetaResolution = 1;
	double rhoResolution = 1;
	// check them
	double rho_max = floor(sqrt(pow(mrows-1,2)+pow(ncols-1,2))+0.5);
    if(thetaResolution>5||rhoResolution>rho_max/10) 
	{
        mexErrMsgIdAndTxt("MyToolbox:voteDCHT:meaningless","Resolution of Theta or Rho is meaningless.");
    }
// --> Prepare data
	
	// compute quantity and step of "theta" & "rho"
	unsigned int thetaQ = (unsigned int)ceil(45/thetaResolution);
	unsigned int rhoQ = (unsigned int)ceil(rho_max/rhoResolution);
	double theta_step = (double)45/thetaQ;
	double rho_step = (double)rho_max/rhoQ;
	// create the output matrix "H" and pass the pointer
    plhs[0] = mxCreateNumericMatrix (2*rhoQ+1,4*thetaQ,mxUINT32_CLASS,mxREAL);
    unsigned int *H = mxGetData(plhs[0]);
	// create the output matrix "theta" and pass the pointer
    plhs[1] = mxCreateDoubleMatrix (1,4*thetaQ,mxREAL);
    double *theta = mxGetPr(plhs[1]);
	// create the output matrix "rho" and pass the pointer
    plhs[2] = mxCreateDoubleMatrix (1,2*rhoQ+1,mxREAL);
    double *rho = mxGetPr(plhs[2]);
	// dynamic memory allocation of "sinMap" and "cosMap"
    double *sinMap = mxMalloc(4*thetaQ*sizeof(double));
    double *cosMap = mxMalloc(4*thetaQ*sizeof(double));
	// compute "theta" & "sinMap" & "cosMap"
	for (unsigned int count = 0; count<thetaQ; count++)
	{
		unsigned int start[4] = {0, thetaQ, 2*thetaQ, 3*thetaQ};
		// compute theta
		theta[start[0]+count] = -90 + theta_step*count;
		theta[start[1]+count] = -45 + theta_step*count;
		theta[start[2]+count] =  0  + theta_step*count;
		theta[start[3]+count] =  45 + theta_step*count;
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
		rho[rhoQ+count] = rho_step*(double)count;
		rho[rhoQ-count] = -rho_step*(double)count;
	}
// --> Main computation

    // call the computational routine
    voteDCHT(img,H,mrows,ncols,rho_step,sinMap,cosMap,thetaQ,rhoQ);
// --> Rest

	// release memory
	mxFree(sinMap);
	mxFree(cosMap);
}
