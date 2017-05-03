/* ===============================================================
 * Copyright (C) 2017, Yirami  
 * Released under BSD License
 * All rights reserved.
 * 
 * Filename: voteSHT.c (compiled into mex)
 * F-Description: vote for Standard Hough transform(SHT)
 *
 *		input 	a MxN matrix "img"
 * 			and 
 *		outputs a PxQ matrix "H"
 *				a 1xR matrix "theta"
 *				a 1xR matrix "rho"
 *
 * 		The calling syntax is:
 *			[H,theta,rho] = voteSHT(img)
 * 
 * Author：Yirami
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
void voteSHT(unsigned int *img, unsigned int *vote, size_t m, size_t n, double rho_step, double *sinMap, double *cosMap, unsigned int thetaQ, unsigned int rhoQ)
{
// --> Vote for each point not zero
	// "j" means cols-1
    for (unsigned int j=0; j<n; j++) 
	{
		// and "i" means rows-1
		for (unsigned int i=0; i<m; i++)
		{
			if (img[j*m+i])
			{
				for (unsigned int num_c=0; num_c<4*thetaQ; num_c++)
				{
					// compute the rho in current point
					double rho_c = floor((j*cosMap[num_c]+i*sinMap[num_c])/rho_step+0.5);
					// vote 
					vote[num_c*(2*rhoQ+1)+(unsigned int)rho_c+rhoQ]++;
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
        mexErrMsgIdAndTxt("MyToolbox:voteSHT:nrhs","One input required.");
    }
    if(nlhs!=3) 
	{
        mexErrMsgIdAndTxt("MyToolbox:voteSHT:nlhs","Three outputs required.");
    }
    // make sure the input argument is type double
    if( !mxIsUint32(prhs[0]) ) 
	{
        mexErrMsgIdAndTxt("MyToolbox:voteSHT:notUint32","Input matrix must be type Uint32.");
    }
    // check that number of dimensions in input argument is 2
    if(mxGetNumberOfDimensions(prhs[0])!=2) 
	{
        mexErrMsgIdAndTxt("MyToolbox:voteSHT:not2D","Input must be 2-D.");
    }
// --> Accept input parameters

    // create a pointer to the real data in the input matrix
    unsigned int *img = mxGetData(prhs[0]);
    // and get dimensions of the input matrix 
	size_t mrows = mxGetM(prhs[0]);
    size_t ncols = mxGetN(prhs[0]);
// --> Parameters can be added to function interface

	// default parameters
	double thetaResolution = 1;
	double rhoResolution = 1;
	// check them
	double rho_max = floor(sqrt(pow(mrows-1,2)+pow(ncols-1,2))+0.5);
    if(thetaResolution>5||rhoResolution>rho_max/10) 
	{
        mexErrMsgIdAndTxt("MyToolbox:voteSHT:meaningless","Resolution of Theta or Rho is meaningless.");
    }
// --> 	Prepare data
	
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
    voteSHT(img,H,mrows,ncols,rho_step,sinMap,cosMap,thetaQ,rhoQ);
// --> Rest

	// release memory
	mxFree(sinMap);
	mxFree(cosMap);
}
