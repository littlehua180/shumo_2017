
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <stdlib.h>

#include "kernel_fun.h"

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )
{
   long i, j, n1, n2;
   double tmp;
   double *K;
  
   if( nrhs == 3) 
   {
      if( !mxIsNumeric(prhs[0]) || !mxIsDouble(prhs[0]) ||
        mxIsEmpty(prhs[0])    || mxIsComplex(prhs[0]) )
        mexErrMsgTxt("Input data must be a real matrix.");

      ker = kernel_id( prhs[1] );
      if( ker == -1 ) 
        mexErrMsgTxt("Improper kernel identifier.");
      
     arg1 = mxGetPr(prhs[2]);

     dataA = mxGetPr(prhs[0]);   
     dataB = dataA;
     dim = mxGetM(prhs[0]);      
     n1 = mxGetN(prhs[0]);       

     plhs[0] = mxCreateDoubleMatrix(n1,n1,mxREAL);
     K = mxGetPr(plhs[0]);

     for( i = 0; i < n1; i++ ) {
        for( j = i; j < n1; j++ ) {
           tmp = kernel( i, j );
           K[i*n1+j] = tmp; 
           K[j*n1+i] = tmp; /* kernel is symetric */
        }
     }
   } 
   else if( nrhs == 4)
   {
      if( !mxIsNumeric(prhs[0]) || !mxIsDouble(prhs[0]) ||
        mxIsEmpty(prhs[0])    || mxIsComplex(prhs[0]) )
        mexErrMsgTxt("Input dataA must be a real matrix.");

      if( !mxIsNumeric(prhs[1]) || !mxIsDouble(prhs[1]) ||
        mxIsEmpty(prhs[1])    || mxIsComplex(prhs[1]) )
        mexErrMsgTxt("Input dataB must be a real matrix.");

      ker = kernel_id( prhs[2] );
      if( ker == -1 ) 
        mexErrMsgTxt("Improper kernel identifier.");

     arg1 = mxGetPr(prhs[3]);

     dataA = mxGetPr(prhs[0]);    
     dataB = mxGetPr(prhs[1]);    
     dim = mxGetM(prhs[0]);       
     n1 = mxGetN(prhs[0]);        
     n2 = mxGetN(prhs[1]);        

     plhs[0] = mxCreateDoubleMatrix(n1,n2,mxREAL);
     K = mxGetPr(plhs[0]);

     for( i = 0; i < n1; i++ ) {
        for( j = 0; j < n2; j++ ) {
           K[j*n1+i] = kernel( i, j );
        }
     }
   }
   else
   {
      mexErrMsgTxt("Wrong number of input arguments.");
   }

   return;
}
