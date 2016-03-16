#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "utilities.h"
#include "power.h"


void PWRComputeError(int n, double *perror, double *newvector, double *vector);

/** free memory of a bag**/
void PWRFreeBag(PowerBag **ppbag)
{
	PowerBag *pbag = *ppbag;

	if (pbag == NULL) goto BACK;

	UTLFree((void**)&pbag->vec0);
	UTLFree((void**)&pbag);

	BACK:
	*ppbag = pbag;
}



int PWRAllocBag(int n, int r, double *covmatrix, PowerBag **ppbag, double scale, double tolerance)
{
	int retcode = 0;
	int i;
	PowerBag *pbag = NULL;
	double *double_array = NULL;

	double *vector, *vector0, *newvector, *q, *qprime, *eigenvalue;

	pbag = (PowerBag *)calloc(1, sizeof(PowerBag));
	if (pbag == NULL) {
		printf("cannot allocate bag\n");
		retcode = NOMEMORY; goto BACK;
	}


	double_array = (double*)calloc(n*r + n*r + n*r + n*n + n*n, sizeof(double));
	if (double_array == NULL) {
		retcode = NOMEMORY; goto BACK;
	}

	/** keep variables used for intensive computations close together in memory for more efficiency **/
	vector0 = &double_array[0];
	vector = &double_array[n*r];
	newvector = &double_array[n*r + n*r];
	qprime = &double_array[n*r + n*r + n*r];
	q = &double_array[n*r + n*r + n*r + n*n];

	/** copy the covariance matrix **/
	for(i = 0; i < n*n; i++)
		q[i] = covmatrix[i];

	eigenvalue = (double*)calloc(n, sizeof(double));

	BACK:
	if (pbag != NULL) {
		/** write bag contents **/
		pbag->n = n;
		pbag->r = r;
		pbag->q = q;
		pbag->qprime = qprime;
		pbag->eigval = eigenvalue;
		pbag->vec = vector;
		pbag->vec0 = vector0;
		pbag->nextvec = newvector;
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
void PWRPCAIteration(int n, double *vector, double *newvector, double *q, double *peigenvalue, double *perror)
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


void PWRPCA(PowerBag *pbag, unsigned int *prseed,
		int (*itercallback)(void *data_itercb, int k, int f, double error), void *data_itercb, int interval_itercb,
		int (*eigvalcallback)(void *data_eigvalcb, int k, int f, double error), void *data_eigvalcb
) {
	int i, j;
	int n, r;
	double *vector, *vector0, *newvector;
	int retcode = 0;
	double tolerance, sp;

	int k, f;
	double error;

	n = pbag->n;
	r = pbag->r;
	vector = pbag->vec;
	vector0 = pbag->vec0;
	newvector = pbag->nextvec;
	tolerance = pbag->tolerance;

	/** initialize first vector to random**/
	for(j = 0; j < n*1; j++){
		vector0[j] = rand_r(prseed)/((double) RAND_MAX);
	}

	/** copy Q into Q'  so that we only deal with Q' and afterwards**/
	for (j = 0; j < n*n; j++)
		pbag->qprime[j] = pbag->q[j];

	for (f = 0; f < r; f++) {
		/** copy f-th column vector0 into vector **/
		for(j = 0; j < n; j++){
			vector[f*n + j] = vector0[f*n + j];
		}
		for(k = 0; ; k++) {

			/* showVector(n, vector);*/
			PWRPCAIteration(n, &vector[f*n], &newvector[f*n], pbag->qprime, &pbag->eigval[f], &error);


			if(error < tolerance) {
				/** finished to compute f-th eigen value **/


				/** Set Q' = Q' - lambda w w^T **/
				for(i = 0; i < n; i++){
					for (j = 0; j < n; j++){
						pbag->qprime[i*n + j] -= pbag->eigval[f]*vector[f*n + i]*vector[f*n + j];
					}
				}

				/** Set w'_0 = w_0 - (w^T w_0) w **/
				if (f < r-1) {
					/** first compute sp = (w^T w_0)**/
					sp = 0.0;
					for (j = 0; j < n; j++) {
						sp += vector[f*n + j] * vector0[f*n + j];
					}
					/** Set w'_0 = w_0 - sp * w **/
					for (j = 0; j < n; j++) {
						vector0[(f+1)*n + j] = vector0[f*n + j] - sp * vector[f*n + j];
					}
				}

				if (eigvalcallback != NULL) {
					retcode = eigvalcallback(data_eigvalcb, k, f, error);
					if (retcode != 0) {
						goto BACK;
					}
				}

				break;
			}

			/**pbag->itercount = k;**/  /** well, in this case we don't really need k **/

			if (k % interval_itercb == 0) {
				if (itercallback != NULL) {
					retcode  = itercallback(data_itercb, k, f, error);
					if (retcode != 0) {
						goto BACK;
					}
				}
			}

		}
	}

	BACK:

	return;
}

