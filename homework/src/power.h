#ifndef POWER_H
#define POWER_H

#define NOMEMORY 100


typedef struct PowerBag {
	int n;
	int r; /** number of eigen values and vectors we want to save in the pca (r == 2 for the homework)**/
	double *q; /** perturbed cov matrix (initially it was called q) used at the begining of a power method iteration**/
	double *qprime; /** Q' cov matrix used in the power method **/
	double *qcopy; /** initial covariance Q (initially it was called matcopy) **/
	double *scratch; /** vector used for the rank 1 perturbation **/
	double *eigenvalue; /** Array of eigen values sorted in decreasing order**/
	double *vector; /** Corresponding matrix of eigen vectors (r x n matrix) **/
	double *vector0; /** Corresponding matrix of eigen vectors at iteration 0 (r x n matrix) **/
	double *newvector; /** Matrix of new eigen vectors at the end of an iteration (r x n matrix)**/
	double scale; /** scale parameter for the rank 1 perturb **/
	double tolerance;
} powerbag;



void PWRFree(void **paddress);
int PWRLoadCov(char *filename, int *pn, double **pmatrix);
int PWRAllocBag(int n, int r, double *covmatrix, powerbag **ppbag, double scale, double tolerance);
void PWRFreeBag(powerbag **ppbag);
void PWRPowerIteration(int n, double *vector, double *newvector, double *q, double *peigenvalue, double *perror);


#endif

