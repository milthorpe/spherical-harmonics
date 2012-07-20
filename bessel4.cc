#include "bessel4.h"

double ftrunc(double x)
{
	if(x >= 0) return floor(x);
	else return ceil(x);
}

double fmax2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? y : x;
}

int imin2(int x, int y)
{
    return (x < y) ? x : y;
}

void J_bessel_HalfInt(double x, int nb, double *b, int *ncalc) {

    int nend, intx, nbmx, i, k, l, m, n, nstart;

    double pold, test;
    double p, alpem, halfx, aa, bb, cc, psave, plast;
    double tover, alp2em, em, en, psavel, sum;

    double sqrtpi = 1.7724538509055159; // sqrt pi

    --b; // Index Shift to Frotran convention
	*ncalc = nb;
	intx = (int) (x);
	for (i = 1; i <= nb; ++i)
	    b[i] = 0.;

	/*===================================================================
	  Branch into  3 cases :
	  1) use 2-term ascending series for small X
	  2) use asymptotic form for large X when NB is not too large
	  3) use recursion otherwise
	  ===================================================================*/

	if (x < rtnsig_BESS) {
	  /* ---------------------------------------------------------------
	     Two-term ascending series for small X.
	     --------------------------------------------------------------- */
	    alpem = 1.5;
	    halfx = (x > enmten_BESS) ? .5 * x :  0.;
	    aa	  =  pow(halfx, 0.5) / (0.5 * sqrtpi) ;
	    bb	  = (x + 1. > 1.)? -halfx * halfx : 0.;
	    b[1] = aa + aa * bb / alpem;

	    if (x != 0. && b[1] == 0.)
	    	*ncalc = 0;
		    /* ----------------------------------------------
		       Calculate higher order functions.
		       ---------------------------------------------- */
		    if (bb == 0.)
		    	tover = (enmten_BESS + enmten_BESS) / x;
		    else
		    	tover = enmten_BESS / bb;

		    cc = halfx;

		    for (n = 2; n <= nb; ++n) {
		    	aa /= alpem;
		    	alpem += 1.;
		    	aa *= cc;
		    	if (aa <= tover * alpem)
		    		aa = 0.;

		    	b[n] = aa + aa * bb / alpem;
		    	if (b[n] == 0. && *ncalc > n)
		    		*ncalc = n - 1;
		    }

	} else if (x > 25. && nb <= intx + 1) {
	    b++; // C convention
	    b[0] = sqrt(2./PI)*sin(x)/sqrt(x); if (nb==1) return;
	    b[1] = sqrt(2./PI)*(sin(x)/x-cos(x))/sqrt(x); if (nb==2) return;
	    for (i=2; i<nb; i++)
	    	b[i] = 2.*(i-.5)/x*b[i-1]-b[i-2];
	    return;

	}
	else { //printf("case 3\n");
    /* rtnsig_BESS <= x && ( x <= 25 || intx+1 < nb ) :
       --------------------------------------------------------
       Use recurrence to generate results.
       First initialize the calculation of P*S.
       -------------------------------------------------------- */
		nbmx = nb - intx;
		n = intx + 1;
		en = (double)(n + n) + 1.;
		plast = 1.;
		p = en / x;
		/* ---------------------------------------------------
		Calculate general significance test.
		--------------------------------------------------- */
		test = ensig_BESS + ensig_BESS;
		if (nbmx >= 3) {
			/* ------------------------------------------------------------
		   Calculate P*S until N = NB-1.  Check for possible overflow.
		   ---------------------------------------------------------- */
			tover = enten_BESS / ensig_BESS;
			nstart = intx + 2;
			nend = nb - 1;
			en = (double) (nstart + nstart) - 1.;
			for (k = nstart; k <= nend; ++k) {
				n = k;
				en += 2.;
				pold = plast;
				plast = p;
				p = en * plast / x - pold;
				if (p > tover) {
				/* -------------------------------------------
			    To avoid overflow, divide P*S by TOVER.
			    Calculate P*S until ABS(P) > 1.
			    -------------------------------------------*/
					tover = enten_BESS;
					p /= tover;
					plast /= tover;
					psave = p;
					psavel = plast;
					nstart = n + 1;
					do {
						++n;
						en += 2.;
						pold = plast;
						plast = p;
						p = en * plast / x - pold;
					} while (p <= 1.);

					bb = en / x;
					/* -----------------------------------------------
				    Calculate backward test and find NCALC,
				    the highest N such that the test is passed.
				    ----------------------------------------------- */
					test = pold * plast * (.5 - .5 / (bb * bb));
					test /= ensig_BESS;
					p = plast * tover;
					--n;
					en -= 2.;
					nend = imin2(nb,n);
					for (l = nstart; l <= nend; ++l) {
						pold = psavel;
						psavel = psave;
						psave = en * psavel / x - pold;
						if (psave * psavel > test) {
							*ncalc = l - 1;
							goto L190;
						}
					}
					*ncalc = nend;
					goto L190;
				}
			}
			n = nend;
			en = (double) (n + n) + 1.;
			/* -----------------------------------------------------
			Calculate special significance test for NBMX > 2.
			-----------------------------------------------------*/
			test = fmax2(test, sqrt(plast * ensig_BESS) * sqrt(p + p));
		}
		/* ------------------------------------------------
		Calculate P*S until significance test passes. */
		do {
			++n;
			en += 2.;
			pold = plast;
			plast = p;
			p = en * plast / x - pold;
		} while (p < test);

L190:
    /*---------------------------------------------------------------
      Initialize the backward recursion and the normalization sum.
      --------------------------------------------------------------- */
		++n;
		en += 2.;
		bb = 0.;
		aa = 1. / p;
		m = n / 2;
		em = (double)m;
		m = (n << 1) - (m << 2);
		/* = 2 n - 4 (n/2)
		   = 0 for even, 2 for odd n */
		if (m == 0)
			sum = 0.;
		else {
			alpem = em - .5;
			alp2em = em + em + .5;
			sum = aa * alpem * alp2em / em;
		}
		nend = n - nb;
		/* if (nend > 0) */
		/* --------------------------------------------------------
        Recur backward via difference equation, calculating
        (but not storing) b[N], until N = NB.
        -------------------------------------------------------- */
		for (l = 1; l <= nend; ++l) {
			--n;
			en -= 2.;
			cc = bb;
			bb = aa;
			aa = en * bb / x - cc;
			m = m ? 0 : 2; /* m = 2 - m failed on gcc4-20041019 */
			if (m != 0) {
				em -= 1.;
				alp2em = em + em + .5;
				if (n == 1)
					break;

				alpem = em - 1. + .5;
				if (alpem == 0.)
					alpem = 1.;
				sum = (sum + aa * alp2em) * alpem / em;
			}
		}
		/*--------------------------------------------------
        Store b[NB].
		--------------------------------------------------*/
		b[n] = aa;
		if (nend >= 0) {
			if (nb <= 1) {
				sum += b[1] * .5;
				goto L250;
			}
			else {
			/*-- nb >= 2 : ---------------------------
		    Calculate and store b[NB-1].
			----------------------------------------*/
				--n;
				en -= 2.;
				b[n] = en * aa / x - bb;
				if (n == 1)
					goto L240;

				m = m ? 0 : 2; /* m = 2 - m failed on gcc4-20041019 */
				if (m != 0) {
					em -= 1.;
					alp2em = em + em + .5;
					alpem = em - .5;
					if (alpem == 0.)
						alpem = 1.;
					sum = (sum + b[n] * alp2em) * alpem / em;
				}
			}
		}

		/* if (n - 2 != 0) */
		/* --------------------------------------------------------
        Calculate via difference equation and store b[N],
        until N = 2.
        -------------------------------------------------------- */
		for (n = n-1; n >= 2; n--) {
			en -= 2.;
			b[n] = en * b[n + 1] / x - b[n + 2];
			m = m ? 0 : 2; /* m = 2 - m failed on gcc4-20041019 */
			if (m != 0) {
				em -= 1.;
				alp2em = em + em + .5;
				alpem = em -  .5;
				if (alpem == 0.)
					alpem = 1.;
				sum = (sum + b[n] * alp2em) * alpem / em;
			}
		}
		/* ---------------------------------------
		Calculate b[1].
		-----------------------------------------*/
		b[1] = 2. * (.5 + 1.) * b[2] / x - b[3];

L240:
		em -= 1.;
		alp2em = em + em + .5;
		if (alp2em == 0.)
			alp2em = 1.;
		sum += b[1] * alp2em;

L250:
		/* ---------------------------------------------------
        Normalize.  Divide all b[N] by sum.
		---------------------------------------------------*/
		/*	    if (.5 + 1. != 1.) poor test */

		sum *= (sqrtpi * pow(.5* x, -.5));

		aa = enmten_BESS;
		if (sum > 1.)
			aa *= sum;
		for (n = 1; n <= nb; ++n) {
			if (fabs(b[n]) < aa)
				b[n] = 0.;
			else
				b[n] /= sum;
		}
	}



}
/*
int main() {
	double x,B[301],fac;
	int nb=301,ncalc;
	for (x=0.00001; x<10000.1; x*=10) {
		J_bessel_HalfInt(x, nb, B, &ncalc);
		fac=sqrt(0.5*PI/x);
		printf("%e %3d - %25.15e %25.15e %25.15e %25.15e\n",x,ncalc,B[0]*fac,B[10]*fac,B[100]*fac,B[300]*fac);
	}
	x=2.248254632739607e+01;
	//J_bessel(&x, &alpha, &nb, B, &ncalc);
	//printf("%d - %.15e\n",ncalc,B[300]*sqrt(PI*.5/x));
	B[300]=besselj(300,x);
	B[299]=besselj(299,x);
	int l;
	for (l=298; l>=0; l--)
		B[l]=(2.*l+3.)*B[l+1]/x-B[l+2];
	for (l=0; l<=300; l++)
		printf("%.15e\n",B[l]);
}*/
