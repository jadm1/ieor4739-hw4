#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>
#include <string.h>
#include <math.h>
#include "utilities.h"
#include "power.h"
#include "mwf.h"

/*
 * todo:
 *
 * add utilities
 * maybe add some utility functions
 *
 * remove the qcopy member from the PowerBag and do the q qcopy stuff in this file (not in power.c)
 * make the small optimization
 *
 * make commented print code able to work if necessary (maybe add a verbosity flag)
 *
 * make code and compilation process portable with windows and mac
 *
 * clean the frankwolfe hw code
 */



int cheap_rank1_perturb(int n, double *scratch, double *matcopy, double *matrix, unsigned int* pseed, double scale);
void show_vector(int n, double *vector);


int master_job_init(int jobnumber, pthread_mutex_t *poutputmutex, void* databag);
int master_job_end(int jobnumber, pthread_mutex_t *poutputmutex, void* databag);
int worker_job(int jobnumber, pthread_mutex_t *poutputmutex, void* databag, int threadID);


int master_job_init(int jobnumber, pthread_mutex_t *poutputmutex, void* databag) {
	/**powerbag *pbag = (powerbag *)databag;**/


	return 0;
}

int master_job_end(int jobnumber, pthread_mutex_t *poutputmutex, void* databag) {
	powerbag *pbag = (powerbag *)databag;
	int j, r;

	r = pbag->r;


	pthread_mutex_lock(poutputmutex);
	for (j = 0; j < r; j++) {
		printf("Job %d: %d-th EigenValue estimate: %g\n", jobnumber, j+1, pbag->eigenvalue[j]);
	}
	/**for (j = 0; j < r; j++) {
		printf("EigenVector #%d: ", j+1);
		PWRshowvector(n, &pbag->vector[j*n]);
	} **//** don't print the eigen vectors they take much room**/
	pthread_mutex_unlock(poutputmutex);

	return 0;
}

int worker_job(int jobnumber, pthread_mutex_t *poutputmutex, void* databag, int threadID) {
	powerbag *pbag = (powerbag *) databag;
	int i, j;
	int n, r;
	int f;
	double *vector, *vector0, *newvector;
	int k;
	int retcode = 0;
	double error, tolerance, sp;
	unsigned int rseed;


	rseed = jobnumber;

	n = pbag->n;
	r = pbag->r;
	vector = pbag->vector;
	vector0 = pbag->vector0;
	newvector = pbag->newvector;
	tolerance = pbag->tolerance;


	/** let's do the perturbation here **/
	/** Q is initialized from qcopy at this line **/
	if((retcode = cheap_rank1_perturb(n, pbag->scratch, pbag->qcopy, pbag->q, &rseed, pbag->scale)))
		goto BACK;


	/** initialize first vector to random**/
	for(j = 0; j < n*1; j++){
		vector0[j] = rand_r(&rseed)/((double) RAND_MAX);
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
			PWRPowerIteration(n, &vector[f*n], &newvector[f*n], pbag->qprime, &pbag->eigenvalue[f], &error);

			if(0 == k%100){
				pthread_mutex_lock(poutputmutex);
				printf("Job %d at iteration %d, %d-th EigenValue: %g, ", jobnumber, k, f+1, pbag->eigenvalue[f]);
				printf("  L1(error) = %.9e\n", error);
				pthread_mutex_unlock(poutputmutex);
			}


			if(error < tolerance) {
				/** finished to compute f-th eigen value **/


				/** Set Q' = Q' - lambda w w^T **/
				for(i = 0; i < n; i++){
					for (j = 0; j < n; j++){
						pbag->qprime[i*n + j] -= pbag->eigenvalue[f]*vector[f*n + i]*vector[f*n + j];
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


				pthread_mutex_lock(poutputmutex);
				printf(" Job %d converged with tol %g! at iter %d. %d-th EigenValue: %g\n", jobnumber, tolerance, k, f+1, pbag->eigenvalue[f]);
				pthread_mutex_unlock(poutputmutex);

				break;
			}

			/**pbag->itercount = k;**/  /** well, in this case we don't really need k **/

			if(0 == k%100) {
				if (MWFWorkerJobCheckInterrupt(threadID) != 0) {
					/** program is interrupted stop working **/
					pthread_mutex_lock(poutputmutex);
					printf(" ID %d interrupting\n", threadID);
					pthread_mutex_unlock(poutputmutex);
					goto BACK;
				}
			}
		}
	}



	BACK:
	return retcode;
}



int main(int argc, char *argv[])
{
	int retcode = 0;
	int j;

	int quantity = 1, numworkers = 1;

	powerbag **ppbag = NULL;
	powerbag *pbag;
	double scale = 1.0;
	double *covmatrix;
	int r, n;
	double tolerance;


	/**
	 * Load parameters
	 * **/

	/**unsigned int rseed = 123;**/

	r = 2; /** default number of factors **/
	tolerance = 1e-6; /** default tolerance parameter**/

	if(argc < 2){
		printf(" usage: rpower filename [-s scale] [-q quantity] [-w workers] [-r num of eigen vals] [-t tolerance]\n");
		retcode = 1; goto BACK;
	}

	for(j = 2; j < argc; j++){
		if (0 == strcmp(argv[j],"-s")){
			j += 1;
			scale = atof(argv[j]);
		}
		else if (0 == strcmp(argv[j],"-q")){
			j += 1;
			quantity = atoi(argv[j]);
		}
		else if (0 == strcmp(argv[j],"-w")){
			j += 1;
			numworkers = atoi(argv[j]);
		}
		else if (0 == strcmp(argv[j], "-r")){
			j += 1;
			r = atoi(argv[j]); /** number of eigen values we want to extract from the PCA on the cov matrix **/
		}
		else if (0 == strcmp(argv[j],"-t")){
			j += 1;
			tolerance = atof(argv[j]);
		}
		else{
			printf("bad option %s\n", argv[j]); retcode = 1; goto BACK;
		}
	}

	if ( numworkers > quantity ){
		numworkers = quantity;
		printf(" --> reset workers to %d\n", numworkers);
	}

	printf("will use scale %g and quantity %d: %d workers, %d eigen values, tolerance: %g\n", scale, quantity, numworkers, r, tolerance);

	/** Load power bags
	 * **/

	ppbag = (powerbag **)calloc(numworkers, sizeof(powerbag *));
	if(!ppbag){
		printf("could not create power bag array\n"); retcode = NOMEMORY; goto BACK;
	}

	retcode = PWRLoadCov(argv[1], &n, &covmatrix); /** read the data once **/
	if (retcode != 0)
		goto BACK;

	for (j = 0; j < numworkers; j++) {
		if ((retcode = PWRAllocBag(n, r, covmatrix, &ppbag[j], scale, tolerance)))
			goto BACK;
	}


	/**
	 * Start the master worker framework
	 */

	MWFMasterThread(quantity, numworkers, (void**)ppbag, master_job_init, master_job_end, worker_job);


	/**
	 * Free power bags
	 */
	for(j = 0; j < numworkers; j++){
		pbag = ppbag[j];
		PWRFreeBag(&pbag);
	}
	UTLFree((void**)&ppbag);

	BACK:
	if (covmatrix != NULL) {
		free(covmatrix); covmatrix = NULL;
	}
	return retcode;
}



int cheap_rank1_perturb(int n, double *scratch, double *matcopy, double *matrix, unsigned int* pseed, double scale)
{
	int retcode = 0, j, i;
	double sum2, invnorm;

	/** first, create a random vector **/
	for(j = 0; j < n; j++)
		scratch[j] = ((double) rand_r(pseed))/((double) RAND_MAX);

	/** next, convert to norm 1 **/
	sum2 = 0;
	for(j = 0; j < n; j++)
		sum2 += scratch[j]*scratch[j];

	invnorm = 1/sqrt(sum2);

	/** rescale **/
	invnorm *= scale;
	for(j = 0; j < n; j++)
		scratch[j] *= invnorm;


	printf("scale for random perturbation: %g\n", scale);

	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			matrix[i*n + j] = scratch[i]*scratch[j] + matcopy[i*n + j];

	return retcode;
}


/** print vector **/
void show_vector(int n, double *vector)
{
	int j;

	for (j = 0; j < n; j++){
		printf(" %g", vector[j]);
	}
	printf("\n");
}




