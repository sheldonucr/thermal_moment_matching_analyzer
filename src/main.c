#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include <iostream>
using namespace std;



int main(int argc, char** argv)
{
        int order = 4; // default order
	int c;
	int optcount = 0;
	char *inputfile;

	if (argc == 1){ // if user type "tmm"
	    cout << "Usage: tmm -i <inputfile> -q <order>\n" << "If order is not assigned, it will be set to default order " << order << ".\n";
            cout << "For example, typing \"tmm -i p4northwood.mf -q 5\" means that we perform the analysis using the input file p4northwood.mf with reduction order 5.\n";
	    return -1;
	}
	while((c = getopt (argc, argv, "hi:q:")) != -1){
	  switch (c)
	    {
	    case 'h': // help
	      cout << "Usage: tmm -i <inputfile> -q <order>\n" << "If order is not assigned, it will be set to default order " << order << ".\n";
              cout << "For example, typing \"tmm -i p4northwood.mf -q 5\" means that we perform the analysis using the input file p4northwood.mf with reduction order 5.\n";
	      optcount ++;
	      return -1;
	    case 'i': // input
	      inputfile = optarg;
	      optcount = optcount + 2;
	      break;
	    case 'q': // order
	      order = atoi(optarg);
	      optcount = optcount +2;
	      break;
	    case '?': // otherwise
              cout << "Type \"tmm\" or \"tmm -h\" for help.\n";
	      optcount ++;
	      return -1;
	    }
	}
	if ((optcount + 1) < argc) { // if user type non-option charactors
	  cout << "Wrong input format.\n";
	  cout << "Usage: tmm -i <inputfile> -q <order>\n" << "If order is not assigned, it will be set to default order " << order << ".\n";
          cout << "For example, typing \"tmm -i p4northwood.mf -q 5\" means that we perform the analysis using the input file p4northwood.mf with reduction order 5.\n";
          return -1;
	}
        if (order < 1 || order > 100){ // if user set inproper input order
          cout << "Bad approximation order, please try again.\n";
	  return -1;
	}

    FILE *infile;
	int i, j;	
	int nscan, matrixRow, matrixCol, ndim, ninput;
	char matrixName, termch;

	const float env = 34.5;
			
	float **G, **C, **B;
	
	/* read in the thermal model in matrices form */
	infile = fopen(inputfile, "r");
	if (infile == NULL){
		printf("%s\n","Matrix file is not found --- try again.");
                cout << "Type \"tmm\" or \"tmm -h\" for help.\n";
		exit(1);		
	}
	else 
	{	
		printf("\n......read matrix file......\n");
	
		while (!feof(infile)){			
 			nscan = fscanf(infile, "%c %i %i%c", &matrixName, &matrixRow, &matrixCol, &termch);
			if (nscan == EOF) break;
			if (nscan != 4 || termch != '\n'){
				nrerror("Improper file format");
			}
						
			switch (matrixName){
			case 'G':
				G = matrix(1, matrixRow, 1, matrixCol);		
				for (i=1; i<=matrixRow; i++) {
					for (j=1; j<=matrixCol; j++) 
						fscanf(infile, "%f ", &G[i][j]);	
				}
				ndim  = matrixRow;
				break;
			case 'C':
				C = matrix(1, matrixRow, 1, matrixCol);
				for (i=1; i<=matrixRow; i++) {
					for (j=1; j<=matrixCol; j++) 
						fscanf(infile, "%f ", &C[i][j]);
				}
				break;
			case 'B':
				B = matrix(1, matrixRow, 1, matrixCol);
				for (i=1; i<=matrixRow; i++) {
					for (j=1; j<=matrixCol; j++) 
						fscanf(infile, "%f ", &B[i][j]);
				}
				ninput = matrixCol;						
				break;
			default:
				nrerror("Unexpected matrix name");
				break;
			}
		}
		printf("\n......finish reading......\n");		
	}
	
    
    /* initialize the tempreture */
	float *initial;
	initial = vector(1, ndim);
	
	for (i=1; i<=ndim; i++) initial[i] = 10.0 + env;	

	
	/* initialize the power mean value and timestep*/
	int nsample = 10;
	float **powertrace, *timestep;
	
	powertrace = matrix(1, ndim, 1, nsample);
	timestep = vector(1, nsample);
	
	for (j=1; j<=nsample; j++) {
		for (i=1; i<=ndim; i++) {		
			powertrace[2][j] = 100.0;
            powertrace[4][j] = 2.5;			
		}
		timestep[j] = 10;		
	}


    /* temperature calculation by tmm */
	float **temperature;
	temperature = matrix(1, ndim, 1, nsample);
	

	tmm(temperature, G, C, env, initial, powertrace, timestep, nsample, ndim, order);
	


	free_matrix(G, 1, ndim, 1, ndim);
	free_matrix(C, 1, ndim, 1, ndim);
	free_matrix(B, 1, ndim, 1, ninput);

	free_matrix(powertrace, 1, ndim, 1, nsample);
	free_matrix(temperature, 1, ndim, 1, nsample);
	free_vector(timestep, 1, nsample);
	
	free_vector(initial, 1, ndim);
	
	return 0;
	
}

