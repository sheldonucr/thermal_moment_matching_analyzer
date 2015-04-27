#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "util.h"
#include <iostream>

using namespace std;
	
/**************************************************************/
/* tmm: thermal moment matching simulation */
/*************************************************************/

void tmm(float **temperature, float **G, float **C,  float env, float *initial, float **powertrace, float *timestep, int nsample, int ndim, int order)
{
	int i, j;
	
	float *pole;
	pole = vector(1, order);
	float simtime = 0.0; // simulation time

    /* assign the node to be mornitored */
	int viewnode = 1;
		
	printf("\n......pole computation......\n");
	
	poleCompute(pole, G, C, viewnode, ndim, order);

	for (i=1; i<=order; i++) {
		printf("pole[%d] = %14.10f\n", i, pole[i]);
	}
	
  
	float **residue, *moment0, *power, temp;
	residue = matrix(1, ndim, 1, order);
	moment0 = vector(1, ndim);
	power = vector(1, ndim);
	
	cout << "\n......temperature computation......\n";
        cout << "\nThermal simulation for node " << viewnode << ".\n";
        cout << "\nnode       time        temperature\n";
	cout << viewnode << "          " << simtime << "          " << initial[viewnode] << "\n";
	
	for(int iter=1; iter<=nsample; iter++) {
	  //	printf("\n......residue computation......\n");

		for (i=1; i<=ndim; i++) {
			power[i] = powertrace[i][iter];
		}
			
	    residueCompute(residue, moment0, pole, G, C, env, initial, power, ndim, order);
	

	    //   printf("\n......temperatur computation......\n");
	
	
		for(i=1; i<=ndim; i++) {		
			temp=0.0;
			for(j=1; j<=order; j++) {
				temp += residue[i][j]*exp(pole[j]*timestep[iter]);
			}
			temperature[i][iter] = temp + moment0[i] + env;
			initial[i] = temperature[i][iter];
			
		      	if (i==viewnode) {
			  simtime = simtime + timestep[iter];
			  cout << viewnode << "          " << simtime << "          " << temperature[i][iter] << "\n";
			  // printf("temperatuer[%d][%d] = %f\n", i, iter, temperature[i][iter]);
			}
		}
	}

	free_matrix(residue, 1, ndim, 1, order);	
	free_vector(pole, 1, order);	
	free_vector(moment0, 1, ndim);
	free_vector(power, 1, ndim);
}
