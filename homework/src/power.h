#ifndef POWER_H
#define POWER_H


typedef struct PowerBag {
	int n;
	int r; /** number of eigen values and vectors we want to save in the pca (r == 2 for the homework)**/
	double *q; /** initial covariance Q **/
	double *qprime; /** Q' cov matrix used in the power method **/
	double *eigval; /** Array of eigen values sorted in decreasing order**/
	double *vec; /** Corresponding matrix of eigen vectors (r x n matrix) **/
	double *vec0; /** Corresponding matrix of eigen vectors at iteration 0 (r x n matrix) **/
	double *nextvec; /** Matrix of new eigen vectors at the end of an iteration (r x n matrix)**/
	double tolerance;
} PowerBag;




int PWRLoadCov(char *filename, int *pn, double **pmatrix);
int PWRAllocBag(int n, int r, double *covmatrix, PowerBag **ppbag, double tolerance);
void PWRFreeBag(PowerBag **ppbag);
void PWRPCAIteration(int n, double *vector, double *newvector, double *q, double *peigenvalue, double *perror);
void PWRPCA(PowerBag *pbag, unsigned int *prseed,
		int (*itercallback)(void *data_itercb, int k, int f, double error), void *data_itercb, int interval_itercb,
		int (*eigvalcallback)(void *data_eigvalcb, int k, int f, double error), void *data_eigvalcb);

#endif

