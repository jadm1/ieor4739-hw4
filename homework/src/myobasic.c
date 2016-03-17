#include <stdio.h>
#include <stdlib.h>
#include "myopt.h"

int myoGetmyoFromFile(myo **ppmyo, char *filename)
{
	int retcode = 0;
	FILE *input = NULL;
	char buffer[100];
	int n, f, j, i;
	myo *pmyo = NULL;

	input = fopen(filename, "r");
	if(!input){
		printf("cannot open file %s\n", filename); retcode = 1; goto BACK;
	}
	printf("reading file %s\n", filename);

	fscanf(input,"%s",buffer);  fscanf(input,"%s",buffer);
	n = atoi(buffer);
	fscanf(input,"%s",buffer);  fscanf(input,"%s",buffer);
	f = atoi(buffer);
	printf("n = %d f = %d\n", n, f);

	if((retcode = myocreatemyo(&pmyo))) goto BACK;

	pmyo->n = n; pmyo->f = f;
	pmyo->lambda = 1.0;

	pmyo->mu = (double *)calloc(n, sizeof(double));
	pmyo->sigma2 = (double *)calloc(n, sizeof(double));
	pmyo->V = (double *)calloc(n*f, sizeof(double));
	pmyo->upper = (double *)calloc(n, sizeof(double));
	pmyo->lower = (double *)calloc(n, sizeof(double));
	pmyo->F = (double *)calloc(f*f, sizeof(double));
	pmyo->gradient = (double *)calloc(n, sizeof(double));
	pmyo->x = (double *)calloc(n, sizeof(double));
	pmyo->Vx = (double *)calloc(f, sizeof(double));
	pmyo->VtF = (double *)calloc(n*f, sizeof(double));

	/** new allocations **/
	pmyo->gradients_and_indices = (grad_sort_struct*)calloc(pmyo->n, sizeof(grad_sort_struct));
	pmyo->gradient_sorted = (double *)calloc(pmyo->n, sizeof(double));
	pmyo->sort_index = (int *)calloc(pmyo->n, sizeof(int));
	pmyo->y_sorted = (double *)calloc(pmyo->n, sizeof(double));
	pmyo->y_best_sorted = (double *)calloc(pmyo->n, sizeof(double));
	pmyo->y = (double *)calloc(pmyo->n, sizeof(double));
	pmyo->Vy = (double *)calloc(pmyo->n, sizeof(double));

	if(!pmyo->mu || !pmyo->sigma2 || !pmyo->V || !pmyo->upper || !pmyo->lower
			|| !pmyo->F || !pmyo->gradient || !pmyo->x || !pmyo->Vx || !pmyo->VtF  ){
		printf("no memory for allocation\n"); retcode = NOMEM; goto BACK;
	}

	fscanf(input,"%s",buffer);  printf("%s\n", buffer);
	for(j = 0; j < n; j++){
		fscanf(input,"%s",buffer);
		pmyo->mu[j] = atof(buffer);
	}

	fscanf(input,"%s",buffer);  printf("%s\n", buffer);
	for(j = 0; j < n; j++){
		fscanf(input,"%s",buffer);
		pmyo->upper[j] = atof(buffer);
	}

	fscanf(input,"%s",buffer);  printf("%s\n", buffer);
	for(j = 0; j < n; j++){
		fscanf(input,"%s",buffer);
		pmyo->sigma2[j] = atof(buffer);
	}


	fscanf(input,"%s",buffer);  printf("%s\n", buffer);
	for(j = 0; j < n; j++){
		for(i = 0; i < f; i++){
			fscanf(input,"%s",buffer);
			pmyo->V[i*n + j] = atof(buffer);
		}
	}

	fscanf(input,"%s",buffer);  printf("%s\n", buffer);
	for(j = 0; j < f; j++){
		for(i = 0; i < f; i++){
			fscanf(input,"%s",buffer);
			pmyo->F[i*f + j] = atof(buffer);
		}
	}


	fclose(input);

	*ppmyo = pmyo;

	BACK:
	printf("done reading with code %d\n", retcode);
	return retcode;
}

void myokillmyo(myo **ppmyo)
{
	myo *pmyo = *ppmyo;

	printf("freeing myo at %p\n", (void *) pmyo);
	free(pmyo->mu);
	free(pmyo->sigma2);
	free(pmyo->V);
	free(pmyo->upper);
	free(pmyo->lower);
	free(pmyo->F);
	free(pmyo->gradient);
	free(pmyo->x);
	free(pmyo->Vx);
	free(pmyo->VtF);

	/** new deallocations **/
	free(pmyo->Vy);
	free(pmyo->y);
	free(pmyo->y_best_sorted);
	free(pmyo->y_sorted);
	free(pmyo->sort_index);
	free(pmyo->gradient_sorted);
	free(pmyo->gradients_and_indices);

	free(pmyo);
	*ppmyo = NULL;
}

int myocreatemyo(myo **ppmyo)
{
	int retcode = 0;
	myo *pmyo;

	pmyo = (myo *)calloc(1, sizeof(myo));
	if(!pmyo){
		printf("no memory for myo\n"); retcode = NOMEM; goto BACK;
	}
	*ppmyo = pmyo;

	printf("created myo at %p\n", (void *) pmyo);

	/*  printf(" Note: each myo takes %d bytes\n", (int) sizeof(myo));
  printf(" n is at %p\n", (void *)&(pmyo->n));
  printf(" f is at %p\n", (void *)&(pmyo->f));
  printf(" sigma2 is at %p\n", (void *)&(pmyo->sigma2));
  printf(" V is at %p\n", (void *)&(pmyo->V));
  printf(" upper is at %p\n", (void *)&(pmyo->upper));*/

	BACK:
	return retcode;
}
