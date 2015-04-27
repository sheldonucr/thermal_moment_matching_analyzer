#ifndef _UTIL_H_
#define _UTIL_H_


void nrerror(char error_text[]);

int *ivector(long nl, long nh);

float *vector(long nl, long nh);

float **matrix(long nrl, long nrh, long ncl, long nch);

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);

void free_ivector(int *v, long nl, long nh);

void free_vector(float *v, long nl, long nh);

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);

void mtm(float **c, float **a, float **b, int row, int col, int rc);

void mtv(float *c, float **a, float *b, int row, int rc);

void printmat(char str[], float **a, int row, int col);

void printvec(char str[], float *a, int row);

void ludcmp(float **a, int n, int *indx, float *d);

void lubksb(float **a, int n, int *indx, float b[]);

void jacobi(float **a, int n, float d[], float **v, int *nrot);

void eigsrt(float d[], float **v, int n);

void balanc(float **a, int n);

void elmhes(float **a, int n);

void hqr(float **a, int n, float wr[], float wi[]);

void poleCompute(float *pole, float **G, float **C, int viewnode, int ndim, int order);

void residueCompute(float **residue, float *moment0, float *pole, float **G, float **C,  float env, float *temperature, float *power, int ndim, int order);

void tmm(float **temperature, float **G, float **C,  float env, float *initial, float **powertrace, float *timestep, int nsample, int ndim, int order);

#endif
