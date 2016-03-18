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


/**
 * B=v0: initial portfolio value
 * v: portfolio value
 * x: position
 * q: quantity of asset
 * w: price of all position in a given asset
 * p: price of 1 unit of a given asset
 *
 * wt = q pt
 * xt = wt / vt
 * vt = sum_n wt = sum_n q pt
 *
 * q = w0/p0 = x0 v0 /p0 =  x0  B/p0
 *
 */


typedef struct Portfolio {
	int n; /** number of assets (whose quantities are non 0) **/
	int t; /** number of periods **/
	double B; /** initial portfolio value **/
	double *p; /** Matrix of prices (unperturbed) **/
	double *q; /** Array of inital fixed quantities **/
	double *delta; /** avg changes **/
	double *sigma; /** std changes **/
	double *pf_values; /** pointer to the output array of portfolio values**/
	double *pf_returns; /** pointer to the output array of portfolio returns**/
} Portfolio;

typedef struct WorkerBag {
	int ID; /** ID (for worker threads) **/
	int num_sim; /** number of simulations **/
	int num_workers; /** number of workers **/
	pthread_mutex_t *poutputmutex;
	Portfolio *pf; /** worker's portfolio non shared variables **/
} WorkerBag;


int load_initial_positions(char* filename, double **px, int* pn, int **pindices, double nonzero_threshold);
int load_prices(char* filename, double **pp, int n, int *indices, int *pt, int max_t);
int compute_avg_changes(double *p, int n, int t, double **pdelta);
int compute_std_changes(double *p, int n, int t, double *delta, double **psigma);

int portfolio_create(Portfolio **ppf, int n, int t, double *x, double *p, double *delta, double *sigma, double B, double *pf_values, double *pf_returns);
void portfolio_delete(Portfolio **ppf);
int portfolio_create_array(int number, Portfolio*** pppf, int n, int t, double *x, double *p, double *delta, double *sigma, double B, double *pf_values, double *pf_returns);
void portfolio_delete_array(int number, Portfolio ***pppf);






int portfolio_simulation(int sim, unsigned int *prseed, pthread_mutex_t *poutputmutex, void* databag, int threadID) {
	int retcode = 0;

	Portfolio *pf = (Portfolio*)databag;
	int n;
	int t;
	double *p, *q, *delta, *sigma;
	int i, j, k;
	double pf_v, pf_v_old; /** portfolio value **/
	double pf_return;

	n = pf->n;
	t = pf->t;
	p = pf->p;
	q = pf->q;
	delta = pf->delta;
	sigma = pf->sigma;

	for (j = 0; j < t; j++) {
		pf_v = 0.0;
		for (i = 0; i < n; i++) {
			pf_v += (p[i*t + j] + (sigma[i]*drawnormal_r(prseed) + delta[i])) * q[i]; /** add current value of the i-th asset + the perturbation **/
		}
		if (j > 0) {
			pf_return += (pf_v - pf_v_old) / pf_v_old;
		}
		pf_v_old = pf_v;
	}

	pf_return /= (t - 1);

	pf->pf_values[sim] = pf_v;
	pf->pf_returns[sim] = pf_return;


	return retcode;
}

void* worker(void * arg) {
	WorkerBag *wbag = (WorkerBag*)arg;
	Portfolio *pf = wbag->pf;
	pthread_mutex_t *poutputmutex = wbag->poutputmutex;
	int ID = wbag->ID;
	int num_sim = wbag->num_sim;
	int num_workers = wbag->num_workers;
	int start, end, num_sim_local;
	int sim;
	unsigned int rseed = ID;

	/** spread the simulations evenly over worker threads **/
	start = ID * (num_sim/num_workers);
	if (ID < num_workers-1)
		end = (ID+1) * (num_sim/num_workers) - 1;
	else
		end = num_sim - 1;
	num_sim_local = end - start + 1;

	pthread_mutex_lock(poutputmutex);
	printf("Worker %d started. %d sims assigned: from %d to %d \n", ID, num_sim_local, start, end);
	pthread_mutex_unlock(poutputmutex);

	for (sim = start; sim <= end; sim++) {
		portfolio_simulation(sim, &rseed, poutputmutex, (void*)pf, ID);
		if (sim % (num_sim_local/10) == 0) {
			printf("W %d: simulation %d, portfolio value: %g, avg daily return: %g %%\n", ID, sim, pf->pf_values[sim], pf->pf_returns[sim]*100);
		}
	}


	return NULL;
}



int main(int argc, char **argv) {
	/**
	 * utility variables
	 */
	int retcode = 0;
	int i, j;
	char *x_filename, *p_filename;
	char *pfv_filename, *pfr_filename;
	FILE *pfv_f, *pfr_f;
	double s;

	/**
	 * Program variables
	 */
	int num_sim, num_workers, max_t;
	WorkerBag *wbag = NULL;
	pthread_t *pthread = NULL;
	pthread_mutex_t outputmutex;
	double *x = NULL;
	int *indices = NULL;
	double *p = NULL;
	double *delta = NULL;
	double *sigma = NULL;
	/**
	 * arrays for results
	 */
	double *pf_values = NULL;
	double *pf_returns = NULL;

	int n; /** number of assets **/
	int t; /** number of periods **/
	Portfolio **ppf = NULL;

	/**
	 * parameters
	 */
	int verbose;
	double B;

	/**
	 * Default parameter values
	 */
	verbose = 0;
	num_sim = 1;
	num_workers = 1;
	max_t = 10000;
	B = 1000000000.0;


	/**
	 * Collect parameters from command line
	 */
	if(argc < 5) {
		printf("usage: %s <portfolio positions file> <prices history file> <portfolio values output file> <portfolio returns output file> [-q simulations number] [-w workers] [-p max periods] [-v verbose]\n", argv[0]);
		retcode = 1; goto BACK;
	}
	for(j = 5; j < argc; j++){
		if (0 == strcmp(argv[j], "-v")){
			j += 0;
			verbose = 1;
		}
		else if (0 == strcmp(argv[j],"-q")){
			j += 1;
			num_sim = atoi(argv[j]);
		}
		else if (0 == strcmp(argv[j],"-w")){
			j += 1;
			num_workers = atoi(argv[j]);
		}
		else if (0 == strcmp(argv[j],"-p")){
			j += 1;
			max_t = atoi(argv[j]);
		}
		else{
			printf("bad option %s\n", argv[j]); retcode = 1; goto BACK;
		}
	}

	if (num_workers > num_sim) {
		num_workers = num_sim;
		printf(" --> reset workers to %d\n", num_workers);
	}

	x_filename = argv[1];
	p_filename = argv[2];
	pfv_filename = argv[3];
	pfr_filename = argv[4];



	printf("loading positions from %s\n", x_filename);
	retcode = load_initial_positions(x_filename, &x, &n, &indices, 1e-7);
	if (retcode != 0) {
		printf("positions could not be loaded !\n"); goto BACK;
	}

	printf("Portfolio assets: %d\n", n);

	if (verbose) {
		printf("x:");
		UTLShowVector(n, x);
		s = 0;
		for (j = 0; j < n; j++)
			s += x[j];
		printf("sum of positions : %g %%\n", s*100.0);
	}

	printf("loading prices from %s\n", p_filename);
	retcode = load_prices(p_filename, &p, n, indices, &t, max_t);
	if (retcode != 0) {
		printf("prices could not be loaded !\n"); goto BACK;
	}

	printf("periods: %d\n", t);

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



	/** allocating memory for result arrays **/
	pf_values = (double*)calloc(num_sim, sizeof(double));
	pf_returns = (double*)calloc(num_sim, sizeof(double));
	if (pf_values == NULL || pf_returns == NULL) {
		retcode = NOMEMORY; goto BACK;
	}



	/**
	 * Creating portfolio array for every worker (no shared memory so we need to make numworkers copies of 1 portfolio)
	 */
	retcode = portfolio_create_array(num_workers, &ppf, n, t, x, p, delta, sigma, B, pf_values, pf_returns);
	if (retcode != 0)
		goto BACK;


	/** check initial portfolio value**/
	s = 0.0;
	for (i = 0; i < n; i++) {
		s += p[i*t + 0] * ppf[0]->q[i];
	}
	printf("Initial portfolio value: %g\n", s);


	/** prepare simple multithreading **/
	wbag = (WorkerBag*)calloc(num_workers, sizeof(WorkerBag));
	if (wbag == NULL) {
		retcode = NOMEMORY; goto BACK;
	}
	pthread = (pthread_t *)calloc(num_workers, sizeof(pthread_t));
	if (pthread == NULL) {
		printf("could not create thread array\n");
		retcode = NOMEMORY; goto BACK;
	}

	pthread_mutex_init(&outputmutex, NULL);

	for(j = 0; j < num_workers; j++) {
		wbag[j].ID = j;
		wbag[j].num_sim = num_sim;
		wbag[j].pf = ppf[j];
		wbag[j].poutputmutex = &outputmutex;
		wbag[j].num_workers = num_workers;
		pthread_mutex_lock(&outputmutex);
		printf("Launching thread for worker %d\n", j);
		pthread_mutex_unlock(&outputmutex);
		pthread_create(&pthread[j], NULL, &worker, (void *)&wbag[j]);
	}

	pthread_mutex_lock(&outputmutex);
	printf("Waiting for threads...\n");
	pthread_mutex_unlock(&outputmutex);

	for(j = 0; j < num_workers; j++) {
		pthread_join(pthread[j], NULL);

		pthread_mutex_lock(&outputmutex);
		printf("Thread %d joined ...\n", j);
		pthread_mutex_unlock(&outputmutex);
	}

	pthread_mutex_destroy(&outputmutex);

	printf("saving results...\n");

	pfv_f = fopen(pfv_filename, "w");
	fprintf(pfv_f, "nsim: %d\n", num_sim);
	for (j = 0; j < num_sim; j++) {
		fprintf(pfv_f, "%g\n", pf_values[j]);
	}
	fclose(pfv_f);

	pfr_f = fopen(pfr_filename, "w");
	fprintf(pfr_f, "nsim: %d\n", num_sim);
	for (j = 0; j < num_sim; j++) {
		fprintf(pfr_f, "%g\n", pf_returns[j]);
	}
	fclose(pfr_f);

	printf("freeing memory ...\n");

	BACK:

	portfolio_delete_array(num_workers, &ppf);

	UTLFree((void**)&pf_values);
	UTLFree((void**)&pf_returns);
	UTLFree((void**)&wbag);
	UTLFree((void**)&pthread);
	UTLFree((void**)&sigma);
	UTLFree((void**)&delta);
	UTLFree((void**)&p);
	UTLFree((void**)&indices);
	UTLFree((void**)&x);

	return retcode;
}


int load_initial_positions(char* filename, double **px, int* pn, int** pindices, double nonzero_threshold) {
	int retcode = 0;
	FILE *f;
	char b[100];
	int i, ii;

	int N; /** total number of assets **/
	int n; /** number of assets in the initial portfolio (whose position is non zero)**/
	double *x = NULL;
	int *indices = NULL;

	f = fopen(filename, "r");
	if (f == NULL) {
		retcode = FILEOPENFAIL; goto BACK;
	}

	fscanf(f, "%s", b);
	fscanf(f, "%s", b);
	N = atoi(b);
	x = (double*)calloc(N, sizeof(double));
	indices = (int*)calloc(N, sizeof(int));
	if (x == NULL || indices == NULL) {
		retcode = NOMEMORY; goto BACK;
	}
	for (i = 0; i < N; i++) {
		fscanf(f, "%s", b);
		x[i] = atof(b);
	}
	fclose(f);

	n = 0;
	for (i = 0; i < N; i++) {
		if (x[i] > nonzero_threshold) {
			indices[n] = i;
			n++;
		}
	}

	/** rearrange the quantities using the indices of nonzero quantities **/
	for (ii = 0; ii < n; ii++) {
		x[ii] = x[indices[ii]];
	}

	if (pindices != NULL)
		*pindices = indices;
	if (pn != NULL)
		*pn = n;
	if (px != NULL)
		*px = x;

	BACK:
	if (retcode != 0) {
		UTLFree((void**)&x);
		UTLFree((void**)&indices);
	}
	return retcode;
}


int load_prices(char* filename, double **pp, int n, int *indices, int *pt, int max_t) {
	int retcode = 0;
	FILE *f;
	char b[100];
	int i, ii, j;

	int N, T;
	int t;
	double *p = NULL;

	f = fopen(filename, "r");
	if (f == NULL) {
		retcode = FILEOPENFAIL; goto BACK;
	}

	fscanf(f, "%s", b);
	fscanf(f, "%s", b);
	N = atoi(b);
	fscanf(f, "%s", b);
	fscanf(f, "%s", b);
	T = atoi(b);
	if (T > max_t)
		t = max_t;
	else
		t = T;

	/** skip the dates line **/
	fscanf(f, "%s", b);
	for (j = 0; j < T; j++)
		fscanf(f, "%s", b);

	p = calloc(n*t, sizeof(double));
	if (p == NULL) {
		retcode = NOMEMORY; goto BACK;
	}
	ii = 0;
	for (i = 0; i < N; i++) {
		fscanf(f, "%s", b);/** jump over asset name **/
		fscanf(f, "%s", b);/** jump over Adj_close: **/
		for (j = 0; j < T; j++) {
			fscanf(f, "%s", b);
			if (ii<n && i == indices[ii] && j < t)
				p[ii*t + j] = atof(b);
		}
		if (i == indices[ii])
			ii++;
	}
	fclose(f);


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



int portfolio_create(Portfolio **ppf, int n, int t, double *x, double *p, double *delta, double *sigma, double B, double *pf_values, double *pf_returns) {
	int retcode = 0;
	int i, j;
	Portfolio *pf = NULL;
	void *memory;

	pf = (Portfolio *)calloc(1, sizeof(Portfolio));
	if (pf == NULL) {
		retcode = NOMEMORY; goto BACK;
	}

	pf->n = n;
	pf->t = t;
	pf->B = B;
	/**
	 * allocate packed memory
	 */
	memory = malloc(n*sizeof(double) +
			n*t*sizeof(double) +
			n*sizeof(double) +
			n*sizeof(double)
	);
	if (memory == NULL) {
		retcode = NOMEMORY; goto BACK;
	}
	pf->q = (double*)memory;
	pf->p = (double*)&pf->q[n];
	pf->delta = (double*)&pf->p[n*t];
	pf->sigma = (double*)&pf->delta[n];


	/**
	 * copy memory
	 */
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
	pf->pf_values = pf_values;
	pf->pf_returns = pf_returns;

	/** compute the initial quantities **/
	for (i = 0; i < n; i++) {
		pf->q[i] = x[i] * (pf->B / pf->p[i*t + 0]);
	}

	if (ppf != NULL)
		*ppf = pf;
	BACK:
	if (retcode != 0) {
		portfolio_delete(&pf);
	}
	return retcode;
}

void portfolio_delete(Portfolio **ppf) {
	Portfolio *pf = *ppf;
	void *memory;
	if (pf == NULL)
		return;

	memory = (void*)pf->q;
	UTLFree(&memory);

	UTLFree((void**)&pf);
	*ppf = NULL;
}


int portfolio_create_array(int number, Portfolio*** pppf, int n, int t, double *x, double *p, double *delta, double *sigma, double B, double *pf_values, double *pf_returns) {
	int retcode = 0;
	int i;
	Portfolio **ppf = NULL;

	ppf = (Portfolio **)calloc(number, sizeof(Portfolio*));
	if (ppf == NULL) {
		retcode = NOMEMORY; goto BACK;
	}

	for (i = 0; i < number; i++) {
		retcode = portfolio_create(&ppf[i], n, t, x, p, delta, sigma, B, pf_values, pf_returns);
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

