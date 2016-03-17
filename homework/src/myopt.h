#ifndef MYO

#define MYO

/** specific structure used to sort gradient components and retrieve the index permutation
  please read myoalgo.c to see what it is used for **/
typedef struct grad_sort_struct {
	double gradient;
	int index;
} grad_sort_struct;

typedef struct myo{
	int n;
	int f;
	int iteration;
	double *mu;
	double *sigma2;
	double *V;
	double *upper;
	double *lower;
	double *F;
	double *gradient;
	double *x;
	double *Vx;
	double *VtF;

	/** variables added to the bag, please read myoalgo.c to understand what they are used for **/
	double lambda;
	double *Vy;
	double *y;
	double *y_best_sorted;
	double *y_sorted;
	int *sort_index;
	double *gradient_sorted;
	grad_sort_struct *gradients_and_indices;
}myo;

#define NOMEM 100
#define DATAERROR1 101

int myoGetmyoFromFile(myo **ppmyo, char *filename);
int myocreatemyo(myo **pmyo);
void myokillmyo(myo **ppmyo);
int myoalgo(myo *pmyo);

#include "utilities.h"

#endif
