extern double *dataA;   /* pointer at the first pattern */
extern double *dataB;   /* pointer at the second pattern */
extern long dim;        /* dimension of patterns */
extern int ker;         /* kernel type (0 - linear, 1 - polynom, 2 - rbf */
extern double *arg1;    /* argument of the kernel */
extern long ker_cnt;    /* number of kernel evaluations */
extern char *kernel_name[]; /* kernel names */


extern double kernel( long a, long b );

int kernel_id( const mxArray *prhs1 );
