#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "utilities.h"
#include "power.h"


void PWRComputeError(int n, double *perror, double *newvector, double *vector);

/** free memory of a bag**/
void PWRFreeBag(powerbag **ppbag)
{
	powerbag *pbag = *ppbag;

	if (pbag == NULL) goto BACK;

	UTLFree((void**)&pbag->vector0);
	UTLFree((void**)&pbag->qcopy);
	UTLFree((void**)&pbag);

	BACK:
	*ppbag = pbag;
}



int PWRAllocBag(int n, int r, double *covmatrix, powerbag **ppbag, double scale, double tolerance)
{
	int retcode = 0;
	int i;
	powerbag *pbag = NULL;
	double *double_array = NULL;

	double *vector, *vector0, *newvector, *q, *qprime, *qcopy, *scratch, *eigenvalue;

	pbag = (powerbag *)calloc(1, sizeof(powerbag));
	if (pbag == NULL) {
		printf("cannot allocate bag\n");
		retcode = NOMEMORY; goto BACK;
	}


	double_array = calloc(n*r + n*r + n*r + n*n + n*n, sizeof(double));
	if (double_array == NULL) {
		retcode = NOMEMORY; goto BACK;
	}

	/** keep variables used for intensive computations close together in memory for more efficiency **/
	vector0 = &double_array[0];
	vector = &double_array[n*r];
	newvector = &double_array[n*r + n*r];
	qprime = &double_array[n*r + n*r + n*r];
	q = &double_array[n*r + n*r + n*r + n*n];

	/** now, allocate an extra matrix and a vector to use in perturbation **/
	/** should really do it in the power retcode since we are putting them in the bag **/
	qcopy = (double *)calloc(n*n, sizeof(double));
	scratch = (double *)calloc(n, sizeof(double));
	if ((qcopy == NULL) || (scratch == NULL)) {
		retcode = NOMEMORY; goto BACK;
	}
	/** and copy the covariance matrix **/
	for(i = 0; i < n*n; i++)
		qcopy[i] = covmatrix[i];

	eigenvalue = (double*)calloc(n, sizeof(double));

	BACK:
	if (pbag != NULL) {
		/** write bag contents **/
		pbag->n = n;
		pbag->r = r;
		pbag->q = q;
		pbag->qprime = qprime;
		pbag->qcopy = qcopy;
		pbag->scratch = scratch;
		pbag->scale = scale;
		pbag->eigenvalue = eigenvalue;
		pbag->vector = vector;
		pbag->vector0 = vector0;
		pbag->newvector = newvector;
		pbag->tolerance = tolerance;
	}
	if (retcode != 0) {
		/** an error occured, cleanup **/
		PWRFreeBag(&pbag);
	}
	*ppbag = pbag;

	return retcode;
}

/** Changed this function to only read the matrix and n from the file
 * This function returns the size of the cov matrix in *pn and the matrix (n*n)
 *
 * **/
int PWRLoadCov(char *filename, int *pn, double **pmatrix)
{
	int retcode = 0, n, j;
	FILE *input = NULL;
	char buffer[100];
	double *matrix = NULL;

	input = fopen(filename, "r");
	if(!input){
		printf("cannot open file %s\n", filename); retcode = 1; goto BACK;
	}
	fscanf(input,"%s", buffer);
	fscanf(input,"%s", buffer);
	n = atoi(buffer);
	printf("n = %d\n", n);


	matrix = (double*)calloc(n*n, sizeof(double));
	if (matrix == NULL) {
		retcode = NOMEMORY;
		goto BACK;
	}

	fscanf(input, "%s", buffer);
	for(j = 0; j < n*n; j++){
		fscanf(input,"%s", buffer);
		matrix[j] = atof(buffer);
	}

	fclose(input);

	BACK:
	if (retcode == 0) {
		printf("read and loaded data for n = %d with code %d\n", n, retcode);
	}
	else {
		UTLFree((void**)&matrix);
		printf("failed to read/load data\n");
	}
	*pn = n;
	*pmatrix = matrix;
	return retcode;
}



/** Compute a power method iteration **/
void PWRPowerIteration(
		int n, double *vector, double *newvector, double *q,
		double *peigenvalue, double *perror)
{
	double norm2 = 0, mult, error;
	int i, j;

	/** w_k+1 = Q * w_k **/
	for(i = 0; i < n; i++){
		newvector[i] = 0;
		for (j = 0; j < n; j++) {
			newvector[i] += vector[j]*q[i*n + j];
		}
	}

	norm2 = 0;
	for(j = 0; j < n; j++)
		norm2 += newvector[j]*newvector[j];

	mult = 1.0/sqrt(norm2);

	for(j = 0; j < n; j++)
		newvector[j] = newvector[j]*mult;

	*peigenvalue = 1.0/mult;

	PWRComputeError(n, &error, newvector, vector);

	/** will need to map newvector into vector if not terminated **/
	for(j = 0; j < n; j++)
		vector[j] = newvector[j];

	*perror = error;
}


void PWRComputeError(int n, double *perror, double *newvector, double *vector)
{
	int j;
	double error;

	error = 0;

	for (j = 0; j < n; j++) {
		error += fabs(newvector[j] - vector[j]);
	}
	error /= n;

	*perror = error;

}



