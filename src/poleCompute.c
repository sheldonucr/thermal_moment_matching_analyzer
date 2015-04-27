#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "util.h"

	
/**************************************************************/
/* poleCompute: compute the poles */
/*************************************************************/


void poleCompute(float *pole, float **G, float **C, int viewnode, int ndim, int order)
{
	int i, j;
	const float etol = 1.0e-15;

	float **R, **RT;
    R = matrix(1, ndim, 1, order);
	RT = matrix(1, order, 1, ndim);

    /* solve G*m_0=b */
    float **a, *b, d;
	int *indx;

	indx = ivector(1, ndim);
	a = matrix(1, ndim, 1, ndim);
	b = vector(1, ndim);

	for (i=1; i<=ndim; i++) {
		for (j=1; j<=ndim; j++) {
			a[i][j] = G[i][j];
		}
		if (i==viewnode) b[i] = 1;
		else b[i] = 0;			   
	}
		
	ludcmp(a, ndim, indx, &d);
	lubksb(a, ndim, indx, b);

    /* normalization */
    float b_sum = 0;
	for (i = 1; i <= ndim; i++){
		b_sum += b[i] * b[i];
	}
	float alpha = sqrt(b_sum);
	for (i = 1; i <= ndim; i++){
		R[i][1] = b[i]/alpha; 
	}


	float csum, h = 0; 

	/* Arnoldi	iteration */		
	for (int iter=2; iter<=order; iter++) {
		for (i=1; i<=ndim; i++) {
        	csum = 0;
			for (j=1; j<=ndim; j++) {
				csum += C[i][j]*R[j][iter-1];
			}			
			b[i] = csum;
		}

		lubksb(a, ndim, indx, b);
		
		for (i=1; i<=ndim; i++) R[i][iter] = b[i];

		for (j=1; j<iter; j++){
			h = 0;
			for (i=1; i<=ndim; i++){
				h += R[i][j] * R[i][iter];
			}
			for (i=1; i<=ndim; i++){
				R[i][iter] -= h * R[i][j];
			}
		}

		b_sum = 0;
		for (i=1; i<=ndim; i++){
			b_sum += R[i][iter] * R[i][iter];
		}
		alpha = sqrt(b_sum);
		if (alpha < etol){
			printf("Break on iteration %d\n", iter);
			order = iter;			
			break;
		}
		else{
			for (i=1; i<=ndim; i++){
				R[i][iter] /= alpha;	
			}
		}
	}
	
	for (j=1; j<=order; j++){
		for (i=1; i<=ndim; i++){			
			RT[j][i] = R[i][j];			
		}
	}
	
    /* model order reduction by congruence transformation */
	float **T, **NG, **NC;
	
	T = matrix(1, order, 1, ndim);
	NG = matrix(1, order, 1, order);
	NC = matrix(1, order, 1, order);
	
	mtm(T, RT, G, order, ndim, ndim);

	mtm(NG, T, R, order, order, ndim);
	

	mtm(T, RT, C, order, ndim, ndim);

	mtm(NC, T, R, order, order, ndim);



	#ifdef DEBUG
    char *mname;
/* 	mname = "R"; */
/* 	printmat(mname, R, ndim, order); */
	mname = "G_r";
	printmat(mname, NG, order, order);
    mname = "C_r";
    printmat(mname, NC, order, order);
	#endif

    #ifdef DEBUG
    printf("\n......eigendecomposition......\n");
    #endif

    /* Eigendcomposition to the reduced model inv(NG)*NC */
	int *indx1;
	float *col, **NGI, **GIC;
	
	indx1 = ivector(1, order);
    col = vector(1, order);
	NGI = matrix(1, order, 1, order);
	GIC = matrix(1, order, 1, order);
	
	ludcmp(NG, order, indx1, &d);
	for (j=1; j<=order; j++) {
		for (i=1; i<=order; i++) col[i]=0.0;
		col[j] = 1.0;
		lubksb(NG, order, indx1, col);
		for (i=1; i<=order; i++) NGI[i][j] = col[i];
	}


	mtm(GIC, NGI, NC, order, order, order);

	#ifdef DEBUG
    mname = "G_r^{-1}";
	printmat(mname, NGI, order, order);
	mname = "G_r^{-1}C_r";
	printmat(mname, GIC, order, order);
	#endif


	float *wr, *wi;
	wr = vector(1, order);
	wi = vector(1, order);
	


	balanc(GIC, order);
	elmhes(GIC, order); 
	hqr(GIC, order, wr, wi); 
    	
    #ifdef DEBUG
    mname = "wr";
	printvec(mname, wr, order);
	#endif
	
/*	jacobi(GIC, order, eig, NGI, &nrot); */

/*	eigsrt(eig, NGI, order); */

	for (i=1; i<=order; i++) wr[i] = fabs(wr[i]);
	
	eigsrt(wr, GIC, order);

    /* For thermal RC circuit, there are only real poles. */

	for (i=1; i<=order; i++) pole[i] = -1.0/wr[i];
	
	#ifdef DEBUG
    mname = "pole";
	printvec(mname, pole, order);
	#endif
	

	free_matrix(R, 1, ndim, 1, order);
	free_matrix(RT, 1, order, 1, ndim);
	free_vector(b, 1, ndim);
	free_matrix(a, 1, ndim, 1, ndim);
	free_matrix(T, 1, order, 1, ndim);
	free_matrix(NG, 1, order, 1, order);
	free_matrix(NC, 1, order, 1, order);
    free_vector(col, 1, order);
	free_matrix(NGI, 1, order, 1, order);
	free_matrix(GIC, 1, order, 1, order);
	free_vector(wr, 1, order);
	free_vector(wi, 1, order);	

	free_ivector(indx, 1, ndim);
	free_ivector(indx1, 1, order);
	
}

