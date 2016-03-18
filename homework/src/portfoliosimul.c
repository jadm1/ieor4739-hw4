/*
 * portfoliosimul.c
 *
 *  Created on: Mar 17, 2016
 *      Author: root
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "utilities.h"




typedef struct Portfolio {
	int n; /** number of assets **/
	int t; /** number of periods **/
	double *p; /** Matrix of prices **/
	double *x; /** Array of positions in each asset**/
	double *delta; /** avg changes **/
	double *sigma; /** std changes **/
	double v; /** portfolio value **/
} Portfolio;


int load_positions(char* filename, double **px, int* pn);
int load_prices(char* filename, double **pp, int* pn, int *pt);
int compute_avg_changes(double *p, int n, int t, double **pdelta);
int compute_std_changes(double *p, int n, int t, double *delta, double **psigma);

int portfolio_create(Portfolio **ppf, int n, int t, double *x, double *p, double *delta, double *sigma);
void portfolio_delete(Portfolio **ppf);
int portfolio_create_array(int number, Portfolio*** pppf, int n, int t, double *x, double *p, double *delta, double *sigma);
void portfolio_delete_array(int number, Portfolio ***pppf);





int master_job_init(int jobnumber, pthread_mutex_t *poutputmutex, void* databag) {


	return 0;
}

int master_job_end(int jobnumber, pthread_mutex_t *poutputmutex, void* databag) {


	return 0;
}

int portfolio_simulation(int jobnumber, pthread_mutex_t *poutputmutex, void* databag, int threadID) {
	int retcode = 0;

	Portfolio *pf = (Portfolio*)databag;
	int n;
	int t;
	unsigned int rseed;

	rseed = jobnumber;

	n = pf->n;
	t = pf->t;




	BACK:
	return retcode;
}



int main(int argc, char **argv) {
	/**
	 * utility variables
	 */
	int retcode = 0;
	int j;
	char *x_filename, *p_filename;
	double s;

	/**
	 * Program variables
	 */
	int quantity, numworkers;
	double *x = NULL;
	double *p = NULL;
	double *delta = NULL;
	double *sigma = NULL;

	int n; /** number of assets **/
	int t; /** number of periods **/
	Portfolio **ppf = NULL;
	Portfolio *pf = NULL;

	/**
	 * parameters
	 */
	int verbose;

	/**
	 * Default parameter values
	 */
	verbose = 0;
	quantity = 1;
	numworkers = 1;



	/**
	 * Collect parameters from command line
	 */
	if(argc < 3) {
		printf("usage: %s <portfolio positions file> <prices history file> [-q quantity] [-w workers] [-v verbose]\n", argv[0]);
		retcode = 1; goto BACK;
	}
	for(j = 3; j < argc; j++){
		if (0 == strcmp(argv[j], "-v")){
			j += 0;
			verbose = 1;
		}
		else if (0 == strcmp(argv[j],"-q")){
			j += 1;
			quantity = atoi(argv[j]);
		}
		else if (0 == strcmp(argv[j],"-w")){
			j += 1;
			numworkers = atoi(argv[j]);
		}
		else{
			printf("bad option %s\n", argv[j]); retcode = 1; goto BACK;
		}
	}

	if (numworkers > quantity) {
		numworkers = quantity;
		printf(" --> reset workers to %d\n", numworkers);
	}

	x_filename = argv[1];
	p_filename = argv[2];

	printf("loading positions from %s\n", x_filename);
	retcode = load_positions(x_filename, &x, &n);
	if (retcode != 0) {
		printf("positions could not be loaded !\n"); goto BACK;
	}

	if (verbose) {
		printf("x:");
		UTLShowVector(n, x);
		s = 0;
		for (j = 0; j < n; j++)
			s += x[j];
		printf("sum of positions : %g %%\n", s*100.0);
	}

	printf("loading prices from %s\n", p_filename);
	retcode = load_prices(p_filename, &p, &n, &t);
	if (retcode != 0) {
		printf("prices could not be loaded !\n"); goto BACK;
	}


	if (verbose) {
		printf("prices loaded\n");
	}

	printf("computing vector of averages of changes...\n");
	retcode = compute_avg_changes(p, n, t, &delta);
	if (retcode != 0)
		goto BACK;
	if (verbose) {
		printf("delta:");
		UTLShowVector(n, delta);
	}

	printf("computing vector std's of changes...\n");
	retcode = compute_std_changes(p, n, t, delta, &sigma);
	if (retcode != 0)
		goto BACK;
	if (verbose) {
		printf("sigma:");
		UTLShowVector(n, sigma);
	}

	/**
	 * Creating portfolio array for every worker (no shared memory so we need to make numworkers copies of 1 portfolio)
	 */
	retcode = portfolio_create_array(numworkers, &ppf, n, t, x, p, delta, sigma);
	if (retcode != 0)
		goto BACK;


	/**
	 * Start the master worker framework
	 */
	MWFMasterThread(quantity, numworkers, (void**)ppf, master_job_init, master_job_end, portfolio_simulation);


	BACK:

	portfolio_delete_array(numworkers, &ppf);

	UTLFree((void**)&sigma);
	UTLFree((void**)&delta);
	UTLFree((void**)&p);
	UTLFree((void**)&x);

	return retcode;
}


int load_positions(char* filename, double **px, int* pn) {
	int retcode = 0;
	FILE *f;
	char b[100];
	int j;

	int n;
	double *x = NULL;

	f = fopen(filename, "r");
	if (f == NULL) {
		retcode = FILEOPENFAIL; goto BACK;
	}

	fscanf(f, "%s", b);
	fscanf(f, "%s", b);
	n = atoi(b);
	x = calloc(n, sizeof(double));
	if (x == NULL) {
		retcode = NOMEMORY; goto BACK;
	}
	for (j = 0; j < n; j++) {
		fscanf(f, "%s", b);
		x[j] = atof(b);
	}
	fclose(f);

	if (pn != NULL)
		*pn = n;
	if (px != NULL)
		*px = x;

	BACK:
	if (retcode != 0) {
		UTLFree((void**)&x);
	}
	return retcode;
}


int load_prices(char* filename, double **pp, int* pn, int *pt) {
	int retcode = 0;
	FILE *f;
	char b[100];
	int i, j;

	int n;
	int t;
	double *p = NULL;



	f = fopen(filename, "r");
	if (f == NULL) {
		retcode = FILEOPENFAIL; goto BACK;
	}

	fscanf(f, "%s", b);
	fscanf(f, "%s", b);
	n = atoi(b);
	fscanf(f, "%s", b);
	fscanf(f, "%s", b);
	t = atoi(b);

	/** skip the dates line **/
	fscanf(f, "%s", b);
	for (j = 0; j < t; j++)
		fscanf(f, "%s", b);

	p = calloc(n*t, sizeof(double));
	if (p == NULL) {
		retcode = NOMEMORY; goto BACK;
	}
	for (i = 0; i < n; i++) {
		fscanf(f, "%s", b);/** jump over asset name **/
		fscanf(f, "%s", b);/** jump over Adj_close: **/
		for (j = 0; j < t; j++) {
			fscanf(f, "%s", b);
			p[i*t + j] = atof(b);
		}
	}
	fclose(f);

	if (pn != NULL)
		*pn = n;
	if (pt != NULL)
		*pt = t;
	if (pp != NULL)
		*pp = p;

	BACK:
	if (retcode != 0) {
		UTLFree((void**)&p);
	}
	return retcode;
}


int compute_avg_changes(double *p, int n, int t, double **pdelta) {
	int retcode = 0;
	int i, k;
	double change;
	double *delta = NULL;

	delta = (double *) calloc(n, sizeof(double));
	if (delta == NULL) {
		retcode = NOMEMORY; goto BACK;
	}

	for (i = 0; i < n; i++) {
		delta[i] = 0;
		for (k = 0; k < (t - 1); k++) {
			change = (p[i * t + k + 1] - p[i * t + k]); /** / p[i * t + k]; (for returns)**/
			delta[i] += change;
		}
		delta[i] /= (t - 1);
	}

	if (pdelta != NULL)
		*pdelta = delta;
	BACK:
	if (retcode != 0) {
		UTLFree((void**)&delta);
	}
	return retcode;
}

int compute_std_changes(double *p, int n, int t, double *delta, double **psigma) {
	int retcode = 0;
	int i, k;
	double change;
	double *sigma = NULL;

	sigma = (double *) calloc(n, sizeof(double));
	if (sigma == NULL) {
		retcode = NOMEMORY; goto BACK;
	}


	for (i = 0; i < n; i++) {
		sigma[i] = 0;
		for (k = 0; k < t - 1; k++) {
			change = (p[i * t + k + 1] - p[i * t + k]);
			sigma[i] += (change - delta[i]) * (change - delta[i]);
		}
		sigma[i] /= (t - 1); /** t - 2 for unbiased estimator**/
		sigma[i] = sqrt(sigma[i]); /** std = sqrt(var) **/
	}


	if (psigma != NULL)
		*psigma = sigma;
	BACK:
	if (retcode != 0) {
		UTLFree((void**)&sigma);
	}
	return retcode;
}



int portfolio_create(Portfolio **ppf, int n, int t, double *x, double *p, double *delta, double *sigma) {
	int retcode = 0;
	int i, j;
	Portfolio *pf = NULL;

	pf = (Portfolio *)calloc(1, sizeof(Portfolio));
	if (pf == NULL) {
		retcode = NOMEMORY; goto BACK;
	}

	pf->n = n;
	pf->t = t;
	/**
	 * allocate memory
	 */
	pf->x = (double*)calloc(n, sizeof(double));
	pf->p = (double*)calloc(n*t, sizeof(double));
	pf->delta = (double*)calloc(n, sizeof(double));
	pf->sigma = (double*)calloc(n, sizeof(double));
	if (pf->x == NULL || pf->p == NULL || pf->delta == NULL || pf->sigma == NULL) {
		retcode = NOMEMORY; goto BACK;
	}

	/**
	 * copy memory
	 */
	for (i = 0; i < n; i++) {
		pf->x[i] = x[i];
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < t; j++) {
			pf->p[i*t + j] = p[i*t + j];
		}
	}
	for (i = 0; i < n; i++) {
		pf->delta[i] = delta[i];
	}
	for (i = 0; i < n; i++) {
		pf->sigma[i] = sigma[i];
	}


	if (ppf != NULL)
		*ppf = pf;
	BACK:
	if (retcode != 0) {
		if (pf != NULL) {
			UTLFree((void**)&pf->x);
			UTLFree((void**)&pf->p);
			UTLFree((void**)&pf->delta);
			UTLFree((void**)&pf->sigma);
		}
		UTLFree((void**)&pf);
	}
	return retcode;
}

void portfolio_delete(Portfolio **ppf) {
	Portfolio *pf = *ppf;
	if (pf == NULL)
		return;

	UTLFree((void**)&pf->x);
	UTLFree((void**)&pf->p);
	UTLFree((void**)&pf->delta);
	UTLFree((void**)&pf->sigma);

	UTLFree((void**)&pf);
	*ppf = NULL;
}


int portfolio_create_array(int number, Portfolio*** pppf, int n, int t, double *x, double *p, double *delta, double *sigma) {
	int retcode = 0;
	int i;
	Portfolio **ppf = NULL;

	ppf = (Portfolio **)calloc(number, sizeof(Portfolio*));
	if (ppf == NULL) {
		retcode = NOMEMORY; goto BACK;
	}

	for (i = 0; i < number; i++) {
		retcode = portfolio_create(&ppf[i], n, t, x, p, delta, sigma);
		if (retcode != 0)
			goto BACK;
	}

	if (pppf != NULL)
		*pppf = ppf;
	BACK:
	if (retcode != 0) {
		portfolio_delete_array(number, pppf);
	}
	return retcode;
}

void portfolio_delete_array(int number, Portfolio ***pppf) {
	int i;
	Portfolio **ppf = *pppf;
	if (ppf == NULL)
		return;

	for (i = 0; i < number; i++)
		portfolio_delete(&ppf[i]);

	UTLFree((void**)&ppf);
	*pppf = NULL;
}

