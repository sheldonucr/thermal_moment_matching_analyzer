#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "util.h"

	
/**************************************************************/
/* residueCompute: compute the residues */
/*************************************************************/


void residueCompute(float **residue, float *moment0, float *pole, float **G, float **C,  float env, float *temperature, float *power, int ndim, int order)
{
	int i, j;

	/* solve moment0 */
    float **a, *b, d;
	int *indx;

	indx = ivector(1, ndim);
	a = matrix(1, ndim, 1, ndim);
	b = vector(1, ndim);

	for (i=1; i<=ndim; i++) {
		for (j=1; j<=ndim; j++) {
			a[i][j] = G[i][j];
		}
		b[i] = power[i];		
	}
		
	ludcmp(a, ndim, indx, &d);
	lubksb(a, ndim, indx, b);

	float **moment, **pm;
    moment = matrix(1, ndim, 1, order);
	pm = matrix(1, order, 1, order);

    /* generate high order moments */

	for (i=1; i<=ndim; i++) {
		moment0[i] = b[i];
		moment[i][1] = moment0[i] - temperature[i] + env;		
	}
	
	float csum; 
	for (int iter=2; iter<=order; iter++) {
		for (i=1; i<=ndim; i++) {
        	csum = 0;
			for (j=1; j<=ndim; j++) {
				csum += C[i][j]*moment[j][iter-1];
			}			
			b[i] = csum;
		}

		lubksb(a, ndim, indx, b);
		
		for (i=1; i<=ndim; i++) moment[i][iter] = -b[i];
	}
	
	for (i=1; i<=order; i++) {
		for (j=1; j<=order; j++) {
			pm[i][j] = -1.0/pow(pole[j],i-1);
		}
	}
	

    /* calculate residues for each node */
	float *t;
	t = vector(1, order);
	int *indx2;
	indx2 = ivector(1, order);
	
	ludcmp(pm, order, indx2, &d);
	
	for (i=1; i<=ndim; i++) {

		for (j=1; j<=order; j++) {
			t[j] = moment[i][j];
		}
		lubksb(pm, order, indx2, t);
		for (j=1; j<=order; j++) {
			residue[i][j] = t[j];
		}
	}
	
	#ifdef DEBUG
	char *mname;	
    mname = "moment0";
	printvec(mname, moment0, ndim);
	#endif

	free_matrix(moment, 1, ndim, 1, order);
	free_vector(b, 1, ndim);
	free_vector(t, 1, order);
	free_matrix(a, 1, ndim, 1, ndim);
	free_matrix(pm, 1, order, 1, order);	

	free_ivector(indx, 1, ndim);
	free_ivector(indx2, 1, order);
	
}
