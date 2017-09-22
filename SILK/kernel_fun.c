#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <string.h>
double *dataA;      /* pointer at the fist patterns */
double *dataB;      /* pointer at the second patterns */
long dim;           /* dimension of patterns */
int ker;            /* kernel id (0 - linear, 1 - polynomial, 
                       2 - rbf, 3 - sigmoid */
double *arg1;       /* kernel argument */
long ker_cnt;       /* number of cernel evaluations */

char *kernel_name[] = {"linear","poly","rbf","sigmoid"};

double dot_prod( long a, long b)
{
   double c = 0;
   long i;
   for( i = 0; i < dim; i++ ) {
      c += *(dataA+(a*dim)+i) * *(dataB+(b*dim)+i);
   }
   return( c );
}

double sub_dot_prod( long a, long b )
{
   double c = 0;
   long i;
   for( i = 0; i < dim; i++ ) {
      c += (*(dataA+(a*dim)+i) - *(dataB+(b*dim)+i))*
           (*(dataA+(a*dim)+i) - *(dataB+(b*dim)+i));
   }
   return( c );
}

int kernel_id( const mxArray *prhs1 )
{
  int num, i, buf_len;
  char *buf;

  if( mxIsChar( prhs1 ) != 1) return( -1 );

  buf_len  = (mxGetM(prhs1) * mxGetN(prhs1)) + 1;
  buf = mxCalloc( buf_len, sizeof( char ));

  mxGetString( prhs1, buf, buf_len );
  
  num = sizeof( kernel_name )/sizeof( char * );

  for( i = 0; i < num; i++ ) {
    if( strcmp( buf, kernel_name[i] )==0 ) return( i );
  }

  return(-1);
}

double kernel( long a, long b )
{
   double c = 0;

   ker_cnt++;

   switch( ker ) {
      case 0:
         c = dot_prod( a, b );
         break;
      case 1:
         c = pow( (dot_prod( a, b) + arg1[1]), arg1[0] );
         break;
      case 2:
         c = exp( -0.5*sub_dot_prod( a, b)/(arg1[0]*arg1[0]) );
         break;
      case 3:     
         c = tanh( arg1[0]*dot_prod( a,b) + arg1[1] );
         break;
      default:
         c = 0;
   }
   return( c );
}

