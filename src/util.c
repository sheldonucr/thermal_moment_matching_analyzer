#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#define NR_END 1
#define FREE_ARG char*

/* #define DEBUG */

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}


float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}


float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}


void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}


void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


void mtm(float **c, float **a, float **b, int row, int col, int rc)
/* matrix times matrix */
{
	int i, j, k;
	
	for (i = 1; i <= row; i++){
		for (j = 1; j <= col; j++){
			c[i][j] = 0;
			for (k = 1; k <= rc; k++){
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}


void mtv(float *c, float **a, float *b, long row, long rc)
/* matrix times vector */
{
	int i, k;
	
	for (i = 1; i <= row; i++){
		c[i] = 0;
		for (k = 1; k <= rc; k++){
			c[i] += a[i][k] * b[k];
		}
	}
}


void printmat(char str[], float **a, int row, int col)
/* print matrix values */
{
	int i, j;
	
	printf("\n");
	
	for (i = 1; i <= row; i++){
		for (j = 1; j <= col; j++){
			printf("%s[%d][%d] = %14.10f\n", str, i, j, a[i][j]);	
		}
	}

	printf("\n");	
}


void printvec(char str[], float *a, int row)
/* print vector values */
{
	int i;
	
	printf("\n");
	
	for (i = 1; i <= row; i++) 
		printf("%s[%d] = %14.10f\n", str, i, a[i]);
	
	printf("\n");
}



/****************************************************************/
/* LU decomposition to solve Ax=b */
/****************************************************************/

# define TINY 1.0e-20;							\
	
void ludcmp(float **a, int n, int *indx, float *d)
/* LU decomposition */
{
	int i, imax, j, k;
	float big, dum, sum, temp;
	float *vv;
	
	vv = vector(1, n);
	*d = 1.0;
	for (i=1; i<=n; i++) {
		big = 0.0;
		for (j=1; j<=n; j++)
			if ((temp=fabs(a[i][j])) > big) big = temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		
		vv[i] = 1.0/big;
	}
	for (j=1; j<=n; j++) {
		for (i=1; i<j; i++) {
			sum = a[i][j];
			for (k=1; k<i; k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i=j; i<=n; i++) {
			sum = a[i][j];
			for (k=1; k<j; k++) 
				sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big = dum;				
				imax = i;
			}
		}
		if (j!=imax) {
			for (k=1; k<=n; k++) {
				dum = a[imax][k];
  				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0)
			a[j][j] = TINY;
		if (j!=n) {
			dum = 1.0/(a[j][j]);
			for (i=j+1; i<=n; i++)
				a[i][j] *= dum;
		}
	}
	free_vector(vv, 1, n);
}


void lubksb(float **a, int n, int *indx, float b[])
/* Forward substitution and backward substitution */
{
	int i, ii=0, ip, j;
	float sum;
	
	for (i=1; i<=n; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii)
			for (j=ii; j<=i-1; j++) 
				sum -= a[i][j]*b[j];
		else if (sum) ii = i;
		b[i] = sum;
	}
	for (i=n; i>=1; i--) {
		sum = b[i];
		for (j=i+1; j<=n; j++)
			sum -= a[i][j]*b[j];
		b[i] = sum/a[i][i];
	}
}

	
/************************************************************/
/* Eigendecomposition by Jacobi method */
/************************************************************/

#define ROTATE(a, i, j, k, l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

void jacobi(float **a, int n, float d[], float **v, int *nrot)
/* Jacobi method */
{
	int i, j, ip, iq;
	float tresh, theta, tau, t, sm, s, h, g, c, *b, *z;
	
	b = vector(1, n);
	z = vector(1, n);
	for (ip=1; ip<=n; ip++) {
		for (iq=1; iq<=n; iq++) v[ip][iq] = 0.0;
		v[ip][ip] = 1.0;
	}
	for (ip=1; ip<=n; ip++) {
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	*nrot = 0;
	for (i=1; i<=50; i++) {
		sm = 0.0;
		for (ip=1; ip<=n-1;ip++) {
			for (iq=ip+1; iq<=n; iq++) 
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_vector(z, 1, n);
			free_vector(b, 1, n);
			return;
		}
		if (i<4) tresh = 0.2*sm/(n*n);
		else tresh = 0.0;
		for (ip=1; ip<=n-1; ip++) {
			for (iq=ip+1; iq<=n; iq++) {
				g = 100.0*fabs(a[ip][iq]);
				if (i>4 && (float)(fabs(d[ip])+g)==(float)fabs(d[ip]) && (float)(fabs(d[iq])+g)==(float)fabs(d[iq]))
					a[ip][iq] = 0.0;
				else if (fabs(a[ip][iq])>tresh) {
					h = d[iq]-d[ip];
					if ((float)(fabs(h)+g)==(float)fabs(h))
						t = (a[ip][iq])/h;
					else {
						theta = 0.5*h/(a[ip][iq]);
						t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta<0.0) t = -t;
					}
					c = 1.0/sqrt(1+t*t);
					s = t*c;
					tau = s/(1.0+c);
					h = t*a[ip][iq];
					z[ip] -=h;
					z[iq] +=h;
					d[ip] -=h;
					d[iq] +=h;
					a[ip][iq] = 0.0;
					for (j=1; j<=ip-1; j++) {
                        ROTATE(a,j,ip,j,iq);
					}
					for (j=ip+1; j<=iq-1; j++) {
						ROTATE(a,ip,j,j,iq);
					}
					for (j=iq+1; j<=n; j++) {
						ROTATE(a,ip,j,iq,j);
					}
					for (j=1; j<=n; j++) {
						ROTATE(v,j,ip,j,iq);
					}
					++(*nrot);
				}
			}
		}
		for (ip=1; ip<=n; ip++) {
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	nrerror("Too many iterations in routine jacobi");
}


/************************************************************/
/* Eigenvalues sorting */
/************************************************************/

void eigsrt(float d[], float **v, int n)
/* sort eigenvalues */
{
	int k,i,j;
	float p;
	
	for (i=1; i<n; i++) {
		p = d[k=i];
		for (j=i+1; j<=n; j++)
			if (d[j]>=p) p = d[k=j];
		if (k!=i) {
			d[k] = d[i];
			d[i] = p;
			for (j=1;j<=n;j++) {
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}


/************************************************************/
/* Balancing a general matrix  */
/************************************************************/
	
#define RADIX 2.0

void balanc(float **a, int n)
{
	int last,i,j;
	float s,r,g,f,c,sqrdx;
	
	sqrdx = RADIX*RADIX;
	last = 0;
	while (last==0) {
		last = 1;
		for (i=1; i<=n; i++) {
			r = c = 0.0;
			for (j=1; j<=n; j++) {
				if (j!=i) {
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
				}
			}
			if (c && r) {
				g = r/RADIX;
				f = 1.0;
				s = c+r;
				while (c<g) {
					f *= RADIX;
					c *= sqrdx;
				}
				g = r*RADIX;
				while (c>g) {
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c+r)/f < 0.95*s) {
					last = 0;
					g = 1.0/f;
					for (j=1; j<=n; j++) a[i][j] *= g;
					for (j=1; j<=n; j++) a[j][i] *= f;
				}
			}
		}
	}
}


#define SWAP(g,h)								\
	{											\
		y = (g);								\
		(g) = (h);								\
		(h) = y;								\
	}											\
		
/************************************************************/
/* reduction to Hessenberg form */
/************************************************************/

void elmhes(float **a, int n)
{
	int m,i,j;
	float x,y;
	
	for (m=2; m<n; m++) {
		x = 0.0;
		i = m;
		for (j=m; j<=n; j++) {
			if (fabs(a[j][m-1])>fabs(x)) {
				x = a[j][m-1];
				i = j;
			}
		}
		if (i!=m) {
			for (j=m-1; j<=n; j++) SWAP(a[i][j],a[m][j]);
			for (j=1; j<=n; j++) SWAP(a[j][i],a[j][m]);
		}
		if (x) {
			for (i=m+1; i<=n; i++) {
				if ((y=a[i][m-1])!=0.0) {
					y /= x;
					a[i][m-1] = y;
					for (j=m; j<=n; j++)
						a[i][j] -= y*a[m][j];
					for (j=1; j<=n; j++) 
						a[j][m] += y*a[j][i];
				}
			}
		}
	}
}


#define SIGN(a, b) ((b)>0.0 ? fabs(a) : -fabs(a))

/************************************************************/
/* QR algorithm for real Hessenberg matrices */
/************************************************************/

void hqr(float **a, int n, float wr[], float wi[])
{
	int nn,m,l,k,j,its,i,nmin;
	float z,y,x,w,v,u,t,s,r,q,p,anorm;
	
	anorm = fabs(a[1][1]);
	for (i=2; i<=n; i++)
		for (j=(i-1); j<=n; j++)
			anorm += fabs(a[i][j]);
	nn = n;
	t = 0.0;
	while (nn>=1) {
		its = 0;
		do {
			for (l=nn; l>=2; l--) {
				s = fabs(a[l-1][l-1])+fabs(a[l][l]);
				if (s==0.0) s = anorm;
				if ((float)(fabs(a[l][l-1])+s)==s) break;
			}
			x = a[nn][nn];
			if (l==nn) {
				wr[nn] = x+t;
				wi[nn--] = 0.0;
			}
			else {
				y = a[nn-1][nn-1];
				w = a[nn][nn-1]*a[nn-1][nn];
				if (l==(nn-1)) {
					p = 0.5*(y-x);
					q = p*p+w;
					z = sqrt(fabs(q));
					x += t;
					if (q>=0.0) {
						z = p + SIGN(z,p);
						wr[nn-1] = wr[nn] = x+z;
						if (z) wr[nn] = x-w/z;
						wi[nn-1] = wi[nn] = 0.0;
					}
					else {
						wr[nn-1] = wr[nn] = x+p;
						wi[nn-1] = -(wi[nn]=z);
					}
					nn -= 2;
				}
				else {
					if (its==30) nrerror("Too many iterations in hqr");
					if (its==10 || its==20) {
						t += x;
						for (i=1; i<=nn; i++) a[i][i] -= x;
						s = fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
						y = x = 0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					for (m=(nn-2); m>=l; m--) {
						z = a[m][m];
						r = x-z;
						s = y-z;
						p = (r*s-w)/a[m+1][m]+a[m][m+1];
						q = a[m+1][m+1]-z-r-s;
						r = a[m+2][m+1];
						s = fabs(p)+fabs(q)+fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m==1) break;
						u = fabs(a[m][m-1])*(fabs(q)+fabs(r));
						v = fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
						if ((float)(u+v)==v) break;
					}
					for (i=m+2; i<=nn; i++) {
						a[i][i-2] = 0.0;
						if (i!=(m+2)) a[i][i-3] = 0.0;
					}
					for (k=m; k<=nn-1; k++) {
						if (k!=m) {
							p = a[k][k-1];
							q = a[k+1][k-1];
							r = 0.0;
							if (k!=(nn-1)) r = a[k+2][k-1];
							if ((x=fabs(p)+fabs(q)+fabs(r))!=0.0) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s=SIGN(sqrt(p*p+q*q+r*r),p))!=0.0) {
							if (k==m) {
								if (l!=m) a[k][k-1] = -a[k][k-1];
							}
							else a[k][k-1] = -s*x;
							p += s;
							x = p/s;
							y = q/s;
							z = r/s;
							q /= p;
							r /= p;
							for (j=k; j<=nn; j++) {
								p = a[k][j]+q*a[k+1][j];
								if (k!=(nn-1)) {
									p += r*a[k+2][j];
									a[k+2][j] -= p*z;
								}
								a[k+1][j] -= p*y;
								a[k][j] -= p*x;
							}
							nmin = nn<k+3 ? nn : k+3;
							for (i=l; i<=nmin; i++) {
								p = x*a[i][k]+y*a[i][k+1];
								if (k != (nn-1)) {
									p += z*a[i][k+2];
									a[i][k+2] -= p*r;
								}
								a[i][k+1] -= p*q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		}
		while (l < nn-1);
	}
}


