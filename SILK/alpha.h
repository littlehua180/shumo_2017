#ifndef _ALPHA_H_
#define _ALPHA_H_

#define  myMAX(a,b)		(a) >= (b) ? (a):(b)
#define  myMIN(a,b)		(a) <= (b) ? (a):(b)

#define MaxReal			1.0e5							/* maximum real value */
#define MinReal			1.0e-15							/* minimum real value */


#define FALSE			0
#define TRUE			1
#define myNULL			0

typedef	struct tagAlpha {								/* component of the self-linked alpha chain */
	double				val;							/* alpha = tao_1 */
	int					id;								/* t \in [0..n-1], where n is the allowed length of the chain */
	struct tagAlpha		*prev, *next;					/* pointers to the previous and next alpha components */
} Alpha;


typedef struct {										/* header of the local alpha chain */
	int					len;							/* # of Alphas currently in the chain */
	Alpha				*first, *last;					/* pointers to the start and end alpha components of the chain */
} Head;

typedef struct tagImgAlpha{								/* alpha chains of the whole image [imh, imw] */
	int					imh, imw;						/* ImSize=imh*imw */
	int					trunc_len;						/* maximal allowed length for the alpha chains: always len <= trunc_len */
	Head				*im;							/* im[ImSize], with im[i] point to one local alpha chain */
} ImgAlpha;

#endif /* _DATA_H_ */

