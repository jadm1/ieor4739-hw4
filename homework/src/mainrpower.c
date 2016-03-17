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



int cheap_rank1_perturb(int n, double *scratch, double *matcopy, double *matrix, unsigned int* pseed, double scale);



int master_job_init(int jobnumber, pthread_mutex_t *poutputmutex, void* databag);
int master_job_end(int jobnumber, pthread_mutex_t *poutputmutex, void* databag);
int pca(int jobnumber, pthread_mutex_t *poutputmutex, void* databag, int threadID);


typedef struct AugPowerBag {
	PowerBag *powerbag;
	double *cov_copy;
	/** cov_copy is a worker local copy of the initial covariance matrix
	 * This additional variable is necessary when we do the perturbation because otherwise the q matrix would stay affected in the powerbag
	 * This is because there is only one PowerBag / worker, not one PowerBag / job
	 * **/
	double *scratch; /** vector used for the rank 1 perturbation **/
	double scale; /** scale parameter for the rank 1 perturb **/
} AugPowerBag;

int main(int argc, char *argv[])
{
	int retcode = 0;
	int i, j;

	int quantity, numworkers;

	AugPowerBag **papbag;
	AugPowerBag *apbag;

	double scale;
	char *cov_filename;
	double *covmatrix = NULL;
	int r, n;
	double tolerance;


	/**
	 * Load parameters
	 * **/

	/** default values **/
	quantity = 1;
	numworkers = 1;
	r = 2; /** default number of factors **/
	tolerance = 1e-6; /** default tolerance parameter**/
	scale = 0.0;

	if(argc < 2){
		printf(" usage: %s filename [-s scale] [-q quantity] [-w workers] [-r num of eigen vals] [-t tolerance]\n", argv[0]);
		retcode = 1; goto BACK;
	}

	cov_filename = argv[1];

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

	/**
	 * Create array of contiguous bags
	 */
	apbag = (AugPowerBag *)calloc(numworkers, sizeof(AugPowerBag) * numworkers);
	/**
	 * Create array of bag pointers (necessary to do that for the MWF function call)
	 */
	papbag = (AugPowerBag **)calloc(numworkers, sizeof(AugPowerBag *));
	if(!papbag){
		printf("could not create bag array\n"); retcode = NOMEMORY; goto BACK;
	}
	for (j = 0; j < numworkers; j++) {
		papbag[j] = &apbag[j]; /** fill it with the corresponding addresses **/
	}


	retcode = PWRLoadCov(cov_filename, &n, &covmatrix); /** read the data once **/
	if (retcode != 0)
		goto BACK;



	/** allocate the power bags **/
	for (j = 0; j < numworkers; j++) {
		apbag = papbag[j];
		if ((retcode = PWRAllocBag(n, r, covmatrix, &apbag->powerbag, tolerance)))
			goto BACK;
	}


	/** specific allocations in the augmented power bag for the rank 1 perturbation **/
	for (j = 0; j < numworkers; j++) {
		apbag = papbag[j];
		apbag->cov_copy = calloc(n*n, sizeof(double));
		if (apbag->cov_copy == NULL) {
			retcode = NOMEMORY; goto BACK;
		}

		/** Save a copy of the covariance matrix in the bag **/
		for (i = 0; i < n*n; i++)
			apbag->cov_copy[i] = covmatrix[i];

		/** now, allocate an extra matrix and a vector to use in perturbation **/
		/** should really do it in the power retcode since we are putting them in the bag **/
		apbag->scratch = (double *)calloc(n, sizeof(double));
		if (apbag->scratch == NULL) {
			retcode = NOMEMORY; goto BACK;
		}

		apbag->scale = scale;
	}



	/**
	 * Start the master worker framework
	 */

	MWFMasterThread(quantity, numworkers, (void**)papbag, master_job_init, master_job_end, pca);


	/**
	 * Free power bags
	 */
	for(j = 0; j < numworkers; j++) {
		apbag = papbag[j];
		PWRFreeBag(&apbag->powerbag);
		UTLFree((void**)&apbag->cov_copy);
	}

	/**
	 * other frees
	 */
	apbag = papbag[0];
	UTLFree((void**)&apbag);
	UTLFree((void**)&papbag);

	UTLFree((void**)&covmatrix);
	
	BACK:
	return retcode;
}






int master_job_init(int jobnumber, pthread_mutex_t *poutputmutex, void* databag) {


	return 0;
}

int master_job_end(int jobnumber, pthread_mutex_t *poutputmutex, void* databag) {
	AugPowerBag *apbag = (AugPowerBag*)databag;
	PowerBag *pbag = apbag->powerbag;
	int j, r;

	r = pbag->r;


	pthread_mutex_lock(poutputmutex);
	for (j = 0; j < r; j++) {
		printf("Job %d: %d-th EigenValue estimate: %g\n", jobnumber, j+1, pbag->eigval[j]);
	}
	/**for (j = 0; j < r; j++) {
		printf("EigenVector #%d: ", j+1);
		PWRshowvector(n, &pbag->vector[j*n]);
	} **//** don't print the eigen vectors they take much room**/
	pthread_mutex_unlock(poutputmutex);

	return 0;
}


typedef struct pcaitercbs {
	pthread_mutex_t *poutputmutex;
	int jobnumber;
	int threadID;
	PowerBag *pbag;
} pcaitercbs;
int pca_itercallback(void *vdata_itercb, int k, int f, double error) {
	pcaitercbs *data_itercb = (pcaitercbs*)vdata_itercb;
	pthread_mutex_t *poutputmutex = data_itercb->poutputmutex;
	int jobnumber = data_itercb->jobnumber;
	int threadID = data_itercb->threadID;
	PowerBag *pbag = data_itercb->pbag;
	int interrupting = 0;


	pthread_mutex_lock(poutputmutex);
	printf("Job %d at iteration %d, %d-th EigenValue: %g, ", jobnumber, k, f+1, pbag->eigval[f]);
	printf("  L1(error) = %.9e\n", error);
	fflush(stdout);
	pthread_mutex_unlock(poutputmutex);


	if (MWFWorkerJobCheckInterrupt(threadID) != 0) {
		/** program is interrupted stop working **/
		pthread_mutex_lock(poutputmutex);
		printf(" ID %d interrupting\n", threadID);
		pthread_mutex_unlock(poutputmutex);
		interrupting = 1;
	}

	return interrupting;
}

typedef struct pcaeigvalcbs {
	pthread_mutex_t *poutputmutex;
	int jobnumber;
	PowerBag *pbag;
} pcaeigvalcbs;
int pca_eigvalcallback(void *vdata_eigvalcb, int k, int f, double error) {
	pcaeigvalcbs *data_eigvalcb = (pcaeigvalcbs*)vdata_eigvalcb;
	pthread_mutex_t *poutputmutex = data_eigvalcb->poutputmutex;
	int jobnumber = data_eigvalcb->jobnumber;
	PowerBag *pbag = data_eigvalcb->pbag;

	pthread_mutex_lock(poutputmutex);
	printf(" Job %d converged with tol %g! at iter %d. %d-th EigenValue: %g\n", jobnumber, pbag->tolerance, k, f+1, pbag->eigval[f]);
	fflush(stdout);
	pthread_mutex_unlock(poutputmutex);

	return 0;
}


int pca(int jobnumber, pthread_mutex_t *poutputmutex, void* databag, int threadID) {
	AugPowerBag *apbag = (AugPowerBag*)databag;
	PowerBag *pbag = apbag->powerbag;
	pcaitercbs data_itercb;
	pcaeigvalcbs data_eigvalcb;

	int n;
	unsigned int rseed;
	int retcode = 0;


	rseed = jobnumber;
	n = pbag->n;

	/** let's do the perturbation here **/
	/** Q is overwritten and reinitialized from a perturbation on cov_copy at this line **/
	if((retcode = cheap_rank1_perturb(n, apbag->scratch, apbag->cov_copy, pbag->q, &rseed, apbag->scale)))
		goto BACK;


	data_itercb.jobnumber = jobnumber;
	data_itercb.threadID = threadID;
	data_itercb.pbag = pbag;
	data_itercb.poutputmutex = poutputmutex;
	data_eigvalcb.jobnumber = jobnumber;
	data_eigvalcb.pbag = pbag;
	data_eigvalcb.poutputmutex = poutputmutex;

	PWRPCA(pbag, &rseed, pca_itercallback, (void*)&data_itercb, 100, pca_eigvalcallback, (void*)&data_eigvalcb);

	BACK:
	/**UTLFree(&qcopy);**/
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







