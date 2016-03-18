/**
 * Code for part A
 * Here we do a PC decomposition using the cov matrix
 */

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>
#include <string.h>
#include <math.h>
#include "utilities.h"
#include "power.h"




int pca_itercallback(void *vpbag, int k, int f, double error) {
	PowerBag *pbag = vpbag;

	printf("Iteration %d, %d-th EigenValue: %g, ", k, f+1, pbag->eigval[f]);
	printf("  L1(error) = %.9e\n", error);
	fflush(stdout);

	return 0;
}

int pca_eigvalcallback(void *vpbag, int k, int f, double error) {
	PowerBag *pbag = vpbag;

	printf("Converged with tol %g! at iter %d. %d-th EigenValue: %g\n", pbag->tolerance, k, f+1, pbag->eigval[f]);
	fflush(stdout);

	return 0;
}


int main(int argc, char *argv[])
{
	int retcode = 0;
	int j;
	int time_ms;

	PowerBag *pbag;

	char *cov_filename;
	char *eigvals_filename;
	char *eigvecs_filename;
	FILE *eigvals_f;
	FILE *eigvecs_f;

	double *covmatrix = NULL;
	int r, n, f;
	double tolerance;
	unsigned int rseed;

	/**
	 * Load parameters
	 * **/

	/** default values **/
	r = 1; /** default number of factors **/
	tolerance = 1e-6; /** default tolerance parameter**/

	if(argc < 4){
		printf(" usage: %s <cov file> <eigvals output file> <eigvecs output file> [-r num of eigen vals] [-t tolerance]\n", argv[0]);
		retcode = 1; goto BACK;
	}

	cov_filename = argv[1];
	eigvals_filename = argv[2];
	eigvecs_filename = argv[3];

	for(j = 4; j < argc; j++){
		if (0 == strcmp(argv[j], "-r")){
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



	printf("will compute %d eigen values, tolerance: %g\n", r, tolerance);


	retcode = PWRLoadCov(cov_filename, &n, &covmatrix); /** read the data once **/
	if (retcode != 0)
		goto BACK;


	/** Load power bag
	 * **/

	/** allocate the power bag **/
	if ((retcode = PWRAllocBag(n, r, covmatrix, &pbag, tolerance)))
		goto BACK;



	/**
	 * PCA
	 */

	rseed = 1;
	time_ms = UTLTicks_ms();
	PWRPCA(pbag, &rseed, pca_itercallback, (void*)pbag, 100, pca_eigvalcallback, (void*)pbag);
	time_ms = UTLTicks_ms() - time_ms;

	for (j = 0; j < r; j++) {
		printf("%d-th EigenValue estimate: %g\n", j+1, pbag->eigval[j]);
	}

	printf("PCA took %d ms\n", time_ms);

	/**
	 * Output results to files
	 */

	eigvals_f = fopen(eigvals_filename, "w");
	fprintf(eigvals_f, "n %d\n", n);
	fprintf(eigvals_f, "r %d\n", r);
	for (f = 0; f < r; f++) {
		fprintf(eigvals_f, "%g\n", pbag->eigval[f]);
	}
	fclose(eigvals_f);

	eigvecs_f = fopen(eigvecs_filename, "w");
	fprintf(eigvecs_f, "n %d\n", n);
	fprintf(eigvecs_f, "r %d\n", r);
	for (f = 0; f < r; f++) {
		for (j = 0; j < n; j++) {
			fprintf(eigvecs_f, "%g ", pbag->vec[f*n + j]);
		}
		fprintf(eigvecs_f, "\n");
	}
	fclose(eigvecs_f);

	/**
	 * Free power bags
	 */
	PWRFreeBag(&pbag);


	/**
	 * other frees
	 */
	UTLFree((void**)&pbag);
	UTLFree((void**)&covmatrix);

	BACK:
	return retcode;
}

