#include <math.h>
#include <stdlib.h>
#include "mex.h"
#include "kernel_fun.h"
#include "alpha.h"


int	 GetLen_alpha_seq(	int r,		int c,
						int imh,	int imw,
						ImgAlpha	*img_alpha)
{
	Head		*pHead;

	pHead					= &(img_alpha->im[c*imh+r]);			/* header of the local alpha chain */

	return( pHead->len );
}

int	 GetVal_alpha_seq(	int r,		int c,	int ind,
						int imh,	int imw,
						ImgAlpha	*img_alpha,
						int*		pId,
						double		decay_rate,
						int			t,
						double*		pVal)
{
	int			i;
	Head		*pHead;
	Alpha		*pCurr;

	pHead					= &(img_alpha->im[c*imh+r]);			/* header of the local alpha chain */
	pCurr					= pHead->first;
	i						= (int)0;
	while (i<ind && pCurr!=(Alpha*)myNULL)
	{
		pCurr				= pCurr->next;
		i					++;
	}

	if (i<ind)
		mexErrMsgTxt("GetVal_alpha_seq():i<ind.\n");
	else
	{
		*pId				= pCurr->id;
		*pVal				= pCurr->val * pow((double)1.0-decay_rate, (double)(t -1 - pCurr->id));
		/* *pVal				= pCurr->val; */
	}

	return ((int)TRUE);
}


int	SetNext_alpha_seq(		int	r,		int c,
						  	int	imh,	int	imw,
							double		val_candid,
							double		decay_rate,
							int			t,
							ImgAlpha	*img_alpha)
{
	Head		*pHead;
	Alpha		*pCurr, *pCand;
	double		val;

	pHead					= &(img_alpha->im[c*imh+r]);			/* header of the local alpha chain */
	pCurr					= pHead->last;

	if (pCurr == (Alpha*)myNULL)
		val					= (double)MaxReal;					/* set to a pseudo-val, if this chain is empty */
	else
		val					= pCurr->val * pow((double)1.0-decay_rate, (double)(t - pCurr->id));

	while (val_candid >= val)
	{
																/* keep upward till val >= val_candid */
		pCurr				= pCurr->prev;
		if (pCurr != (Alpha*)myNULL)
			val				= pCurr->val * pow((double)1.0-decay_rate, (double)(t - pCurr->id));
		else
			break;
	}

	if (pHead->len == img_alpha->trunc_len && pCurr!=(Alpha*)myNULL && pCurr->next ==(Alpha*)myNULL)
		return ((int)FALSE);

	if ((pCand				= (Alpha*) malloc(sizeof(Alpha))) == myNULL)
		mexErrMsgTxt("SetNext_alpha_seq():couldn't allocate memory for (Alpha*)pCand.\n");
	pCand->id				= t;								/* 0-based, since t is also 0-based */
	pCand->val				= val_candid;
	pCand->prev				= pCurr;
	pHead->len				++;
	if (pCurr != (Alpha*)myNULL)
	{
		pCand->next			= pCurr->next;
		pCurr->next			= pCand;
	}
	else														/* pCand is the new first */
	{
		pCand->next			= pHead->first;
		pHead->first		= pCand;
	}
	if (pCand->next !=  (Alpha*)myNULL)
		pCand->next->prev	= pCand;

	if (pCand->next == (Alpha*)myNULL)
																/* pCand is the new last */
		pHead->last			= pCand;

	if (pHead->len > img_alpha->trunc_len)
	{
		pCurr				= pHead->last;
		pCurr->prev->next	= (Alpha*)myNULL;
		pHead->last			= pCurr->prev;
		pHead->len			--;
		free(pCurr);
	}
	
	return ((int)TRUE);
}


int	Init_alpha_seq(			int		imh,	int		imw,
							int			trunc_len,	
							ImgAlpha	*img_alpha)

{
	int r, c;

	img_alpha->imh				= imh;
	img_alpha->imw				= imw;
	img_alpha->trunc_len		= trunc_len;
	if ((img_alpha->im			= (Head*) malloc(imh*imw*sizeof(Head))) == myNULL)
		mexErrMsgTxt("Init_alpha_seq():couldn't allocate memory for img_alpha->im.\n");

	for (r=0;r<imh;r++)
	for (c=0;c<imw;c++)
	{														/* init the empty headers */
		img_alpha->im[c*imh+r].len	= (int) 0;
		img_alpha->im[c*imh+r].first= (Alpha*) myNULL;
		img_alpha->im[c*imh+r].last	= (Alpha*) myNULL;
	}

	return ((int)TRUE);
}


int	Copy2alpha_seq(	int		imh,	int		imw,
					int		bWant_alpha_seq,
					ImgAlpha	*img_alpha, 
					double		*alpha_seq)
{
	int			r, c, t, i;
	Head		*pHead;
	Alpha		*pCurr;

	if (bWant_alpha_seq==(int)0 )
		return ((int)TRUE);

	for (r=0;r<imh;r++)
	for (c=0;c<imw;c++)
	{														/* for each local alpha chain */
		pHead					= &(img_alpha->im[c*imh+r]);
		t						= (int)0;
		pCurr					= pHead->first;
		while (pCurr!=(Alpha*)myNULL)
		{
			if (bWant_alpha_seq==(int)1)
				alpha_seq  [t*imh*imw  +  c*imh+r]
								= pCurr->val;

			pCurr				= pCurr->next;
			t					++;
		}
		for(i=t;i<img_alpha->trunc_len;i++)
		{
			if (bWant_alpha_seq==(int)1)
				alpha_seq  [i*imh*imw  +  c*imh+r]
								= (double)0.0;

		}
	}

	return ((int)TRUE);
}


int	Del_alpha_seq(	int		imh,	int		imw,
					ImgAlpha	*img_alpha)

{
	int			r, c;
	Head		*pHead;
	Alpha		*pCurr, *pCand;

	for (r=0;r<imh;r++)
	for (c=0;c<imw;c++)
	{														/* for each header, delete its local chain */
		pHead					= &(img_alpha->im[c*imh+r]);
		pCurr					= pHead->first;
		while (pCurr != (Alpha*)myNULL)
		{
			pCand				= pCurr;
			pCurr				= pCurr->next;
			free(pCand);
		}
	}
	free(img_alpha->im);

	return ((int)TRUE);
}


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{
	int			imh, imw, T, nTrImgs, nFrames, id, len;
	int			bWant_alpha_seq;
	int			i, j, c, r, t, i_tmp, trunc_len;
	double		radius, C_train, C_tst_large, C_tst_small, threshold, decay_rate;
	double		xt_len2, y_t, ft_xt, l_t, alpha_i;
	double		*trimgs, *imseq, *imseq_seg, *alpha_seq, *imseq_loss;
	double		*x_t, *x_i;
	mxArray		*pr, *subs_xt, *subs_xi, *subs, *subs_imseq_loss, *subs_imseq_seg; 
	ImgAlpha	*img_alpha;

	if (nrhs != 3) { 
		mexErrMsgTxt("Three input arguments required."); 
	} 
	if (nlhs != 0) {
		mexErrMsgTxt("Zero output arguments required."); 
	} 

	pr			= (mxArray *) mexGetVariablePtr("global", "imh");
	imh			= (int) mxGetScalar(pr);
	pr			= (mxArray *) mexGetVariablePtr("global", "imw");
	imw			= (int) mxGetScalar(pr);
	pr			= (mxArray *) mexGetVariablePtr("global", "dim");
	dim			= (int) mxGetScalar(pr);
	pr			= (mxArray *) mexGetVariablePtr("global", "nTrImgs");
	nTrImgs		= (int) mxGetScalar(pr);
	pr			= (mxArray *) mexGetVariablePtr("global", "nFrames");
	nFrames		= (int) mxGetScalar(pr);
	subs		= (mxArray *) mexGetVariablePtr("global", "trimgs");
	trimgs		= mxGetPr(subs);
	subs		= (mxArray *) mexGetVariablePtr("global", "imseq");
	imseq		= mxGetPr(subs);
	pr			= (mxArray *) mexGetVariablePtr("global", "radius");
	radius		= mxGetScalar(pr);
	pr			= (mxArray *) mexGetVariablePtr("global", "C_train");
	C_train		= mxGetScalar(pr);
	pr			= (mxArray *) mexGetVariablePtr("global", "C_tst_large");
	C_tst_large	= mxGetScalar(pr);
	pr			= (mxArray *) mexGetVariablePtr("global", "C_tst_small");
	C_tst_small	= mxGetScalar(pr);
	pr			= (mxArray *) mexGetVariablePtr("global", "threshold");
	threshold	= mxGetScalar(pr);
	pr			= (mxArray *) mexGetVariablePtr("global", "trunc_len");
	trunc_len	= (int) mxGetScalar(pr);
	pr			= (mxArray *) mexGetVariablePtr("global", "tau");
	decay_rate	= mxGetScalar(pr);
	subs		= (mxArray *) mexGetVariablePtr("global", "alpha_seq");
	alpha_seq	= mxGetPr(subs);
	subs_imseq_loss		= (mxArray *) mexGetVariablePtr("global", "imseq_loss");
	imseq_loss	= mxGetPr(subs_imseq_loss);
	subs_imseq_seg		= (mxArray *) mexGetVariablePtr("global", "imseq_seg");
	imseq_seg	= mxGetPr(subs_imseq_seg);


	ker					= kernel_id( prhs[0] );
	if( ker == -1 ) 
		mexErrMsgTxt("Improper kernel identifier.");
	arg1				= mxGetPr(prhs[1]);
	bWant_alpha_seq		= (int)mxGetScalar(prhs[2]);

	subs_xt				= mxCreateDoubleMatrix(dim, 1,   mxREAL);
	x_t					= mxGetPr( subs_xt );
	subs_xi				= mxCreateDoubleMatrix(dim, 1,   mxREAL);
	x_i					= mxGetPr( subs_xi );

	if ((img_alpha				= (ImgAlpha*) malloc(sizeof(ImgAlpha))) == myNULL)
		mexErrMsgTxt("SOLK_1svm_alphaSeq():couldn't allocate memory for (ImgAlpha*)img_alpha.\n");
	Init_alpha_seq(imh, imw, trunc_len, img_alpha);

	T							= nTrImgs;
	for (t=0;t<T;t++)
	{
		for (r=0;r<imh;r++)
		for (c=0;c<imw;c++)
		{
			for (j=0;j<dim;j++)
				x_t[j]			= trimgs[t*(imh*imw*dim)+ j*(imh*imw) +c*imh+r];
			xt_len2				= (double)1.0;
			len					= GetLen_alpha_seq(r,c, imh, imw, img_alpha);
			ft_xt				= (double)0.0;
			for (i=0;i<len;i++)
			{
				GetVal_alpha_seq(r,c,i, imh, imw, img_alpha, &id, decay_rate, 	t, &alpha_i);
				for (j=0;j<dim;j++)
					x_i[j]		= trimgs[id*(imh*imw*dim)+ j*(imh*imw) +c*imh+r];
				dataA			= x_t;
				dataB			= x_i;
				ft_xt			+=	alpha_i * kernel(0, 0);
			}

			l_t					= myMAX(radius-(1-decay_rate)*ft_xt, (double)0.0);

			SetNext_alpha_seq(      r,c,  imh, imw, myMIN( l_t/xt_len2, C_train*(1-decay_rate) ),
									decay_rate, 	t, img_alpha	);

		}

	}

	T							= nFrames;
	for (t=0;t<T;t++)
	{

		for (r=0;r<imh;r++)
		for (c=0;c<imw;c++)
		{
			for (j=0;j<dim;j++)
				x_t[j]			= imseq[t*(imh*imw*dim)+ j*(imh*imw) +c*imh+r];
			xt_len2				= (double)1.0;

			len					= GetLen_alpha_seq(r,c, imh, imw, img_alpha);
			ft_xt				= (double)0.0;
			for (i=0;i<len;i++)
			{
				GetVal_alpha_seq(r,c,i, imh, imw, img_alpha, &id, decay_rate, 	t, &alpha_i);

				if (id<nTrImgs)
				{													/* x_i is in the training imgs */
					i_tmp		= id;
					for (j=0;j<dim;j++)
						x_i[j]	= trimgs[i_tmp*(imh*imw*dim)+ j*(imh*imw) +c*imh+r];
				}
				else
				{													/* x_i is in the running imgs */
					i_tmp		= id - nTrImgs;
					for (j=0;j<dim;j++)
						x_i[j]	= imseq[i_tmp*(imh*imw*dim)+ j*(imh*imw) +c*imh+r];
				}

				dataA			= x_t;
				dataB			= x_i;
				ft_xt			+= alpha_i * kernel(0, 0);


			}
			l_t					= myMAX(radius-(1-decay_rate)*ft_xt, (double)0.0);

			if (l_t >= threshold)
			{
				y_t				= (double)1.0;
				SetNext_alpha_seq(      r,c,  imh, imw, myMIN( l_t/xt_len2, C_tst_small*(1-decay_rate) ),
										decay_rate, 	t+nTrImgs, img_alpha			);
			}
			else
			{
				y_t				= (double)0.0;
				SetNext_alpha_seq(      r,c,  imh, imw, myMIN( l_t/xt_len2, C_tst_large*(1-decay_rate) ),
										decay_rate,  t+nTrImgs, img_alpha			);
			}
			imseq_seg[t*imh*imw  +  c*imh+r]
								= y_t;
			imseq_loss[t*imh*imw  +  c*imh+r]
								= l_t;

		}	


	}	

	Copy2alpha_seq(imh, imw, bWant_alpha_seq, img_alpha, alpha_seq);
	Del_alpha_seq( imh, imw, img_alpha);

    mexPutVariable("global", imseq_loss, subs_imseq_loss);
    mexPutVariable("global", imseq_seg, subs_imseq_seg);

	mxDestroyArray(subs_xi);
	mxDestroyArray(subs_xt);

} 
