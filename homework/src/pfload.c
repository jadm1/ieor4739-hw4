
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utilities.h"
#include "pf.h"



/** local function **/
void portfolio_compute_qty_from_alloc(int n, int t, double B, double *q, double *x, double *p);


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
	portfolio_compute_qty_from_alloc(n, t, pf->B, pf->q, x, pf->p);
	
	if (ppf != NULL)
		*ppf = pf;
	BACK:
	if (retcode != 0) {
		portfolio_delete(&pf);
	}
	return retcode;
}


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
void portfolio_compute_qty_from_alloc(int n, int t, double B, double *q, double *x, double *p) {
	int i;
	
	for (i = 0; i < n; i++) {
		q[i] = x[i] * (B / p[i*t + 0]);
	}
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


