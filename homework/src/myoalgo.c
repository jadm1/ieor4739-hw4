#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myopt.h"


int myofind_feasible(myo *pmyo);
int myo_iteration(myo *pmyo);
int myo_getgradient(myo *pmyo);
int myo_step(myo *pmyo);
void myo_showx(myo *pmyo, int start, int end);
void myoVtimesy(myo *pmyo, double *y);
int myoprepare(myo *pmyo);
int compare_grad_components(const grad_sort_struct* a, const grad_sort_struct* b);

/**#define LOUDFEASIBLE**/

int myoalgo(myo *pmyo)
{
	int retcode = 0;
	int max_its = pmyo->max_iter;

	if((retcode = myoprepare(pmyo))) goto BACK;

	if((retcode = myofind_feasible(pmyo))) goto BACK;

	pmyo->iteration = 0;
	while (1) {
		if ((pmyo->iteration >= 2) && (fabs(pmyo->cost - pmyo->lastcost) < pmyo->mingap)) {
			break;
		}
		if (pmyo->iteration >= max_its) {
			break;
		}

		if((retcode = myo_iteration(pmyo))) goto BACK;
		if (pmyo->verbose)
			myo_showx(pmyo, 0, pmyo->n-1);
		fflush(stdout);

		pmyo->iteration++;
	}

	BACK:
	printf("\nmyoalgo returning with code %d\n\n", retcode);
	return retcode;
}


int myoprepare(myo *pmyo)
{
	int retcode = 0;
	int i, j, k, f = pmyo->f, n = pmyo->n;
	double sum, *V = pmyo->V, *F = pmyo->F;
	/** first compute VtF, which is an nxf matrix **/

	for(j = 0; j < n; j++){
		for(i = 0; i < f; i++){
			/* compute the (j,i) entry of VtF */
			sum = 0;
			for (k = 0; k < f; k++){
				sum  += V[k*n + j]*F[k*f + i];
			}
			pmyo->VtF[j*f + i] = sum;
		}
	}
	printf("myoprepare returns %d\n", retcode);
	return retcode;
}

int myofind_feasible(myo *pmyo)
{
	int retcode = 0;
	int j;
	double sum, *x = pmyo->x, delta;

	sum = 0;
	for(j = 0; j < pmyo->n; j++){
		x[j] = pmyo->lower[j];  /** we assume lower <= upper, but should
				      check **/
#ifdef LOUDFEASIBLE
		printf("  -> x[%d] initialized at %g\n", j, x[j]);
#endif

		sum += x[j];
	}
	if(sum > 1.0){
		printf(" error: sum of asset lower bounds equals %g > 1.0\n", sum);
		retcode = DATAERROR1; goto BACK;
	}

	for(j = 0; (j < pmyo->n) && (sum < 1.0) ; j++){
		if(sum < 1.0){
			delta = (1.0 - sum < pmyo->upper[j] - pmyo->lower[j]) ?
					1.0 - sum : pmyo->upper[j] - pmyo->lower[j];
			sum += delta;
			x[j] += delta;
		}
#ifdef LOUDFEASIBLE
		printf("  -> x[%d] increased to %g\n", j,x[j]);
#endif
	}
	printf("find_feasible done at j = %d\n", j);

	BACK:
	return retcode;
}

int myo_iteration(myo *pmyo)
{
	int retcode = 0;

	if (pmyo->verbose)
		printf("\niteration %d\n", pmyo->iteration);

	if( (retcode = myo_step(pmyo))) goto BACK;

	BACK:
	return retcode;
}

int myo_getgradient(myo *pmyo)
{
	int retcode = 0, i, j, n = pmyo->n, f = pmyo->f;
	double *Vx = pmyo->Vx, *VtF = pmyo->VtF, sum;

	if (pmyo->verbose)
		printf(" computing gradient at iteration %d\n", pmyo->iteration);
	myoVtimesy(pmyo, pmyo->x);
	for(j = 0; j < n; j++){
		/** first initialize the gradient **/

		pmyo->gradient[j] = -pmyo->mu[j] + 2*pmyo->lambda*pmyo->sigma2[j]*pmyo->x[j];

		/** computes jth entry of VtFVx */
		sum = 0;
		for(i = 0; i < f; i++){
			sum += VtF[j*f + i]*Vx[i];
		}
		pmyo->gradient[j] += 2*pmyo->lambda*sum;
	}

	return retcode;
}


int myo_step(myo *pmyo)
{
	/**  **/

	/** The following variables are variables of the myo structure
		So these variables are the important variables of the problem
		Most of them involve allocations, which is why we are defining them before this function to avoid unnecessary mallocs and frees
		Here, we simply define local pointer copies of the bag variables to have a clear code. Otherwise pmyo-> would be everywhere
	 **/
	int n = pmyo->n; /** Number of assets **/
	int f = pmyo->f; /** Number of factors **/

	double* x = pmyo->x; /** Current position whose cost we are trying to optimize **/
	double* y = pmyo->y; /** Feasible direction that we will determine by solving an intermediate linear problem **/
	double* V = pmyo->V; /** Matrix V **/
	double* Vy = pmyo->Vy; /** Product V times y **/
	double* VtF = pmyo->VtF; /** Product V^T times F **/
	double* sigma2 = pmyo->sigma2; /** Diagonal of sigma^2 of the residual matrix D = Q - V^T F V **/
	double* upper = pmyo->upper; /** Upper bounds of x **/
	double* lower = pmyo->lower; /** Lower bound of x **/
	double* gradient = pmyo->gradient; /** Gradient of the cost function **/
	double lambda = pmyo->lambda; /** lambda ( = 1.0 in our problem) **/

	grad_sort_struct* gradients_and_indices = pmyo->gradients_and_indices;
	/** This is an array of n elements with structure grad_sort_struct 
	Each element of that structure is made of a double and an integer.
	The double represents a component of the gradient
	The integer represents the index of that component.
	Doing this will turn out to be useful when calling the quick sort (qsort) function on this array,
	because even if only the doubles will be used to sort the elements, the indices
	will be sorted at the same time as the components.
	This will be useful to get the permutation of indices that was done to sort the components.
	 **/

	double* gradient_sorted = pmyo->gradient_sorted; /** array of the n gradient components sorted in a decreasing order **/
	int* sort_index = pmyo->sort_index; /** corresponding array of the n indices of the permutation done to get the sorted gradient comps **/
	/** Remark : we will have: gradient[sort_index[i]] == gradient_sorted[i] **/

	double* y_sorted = pmyo->y_sorted; /** This is the y used when solving the auxiliary LP (its components correspond to the sorted permutation sort_index) **/
	double* y_best_sorted = pmyo->y_best_sorted; /** Used to store the best y_sorted when enumerating feasible y's in the auxiliary LP **/



	/** These variables are less important and typically hold temporary data for algorithmic purposes **/
	int retcode = 0; /** return code **/
	int i, j; /** counters **/
	int m; /** counter used for the enumeration to determine the best feasible y **/
	double sum; /** temporary value used in calculations (matrix products, etc...) **/
	int y_feasible; /** boolean variable used to check y's feasibility **/
	int num_feasible_y; /** number of feasible y's found after enumeration **/
	double cost; /** cost of the linear problem used to determine y (the direction) (cost == gk^T y) **/
	double bestcost; /** best cost used to compare costs and store the best one found **/
	double s; /** step size **/
	double a, b; /** later we will see that s satisfies a linear equation : as + b = 0, so these variables will be used to store the coefficients **/


	/** compute the gradient **/
	if( (retcode = myo_getgradient(pmyo))) goto BACK;

	if (pmyo->verbose)
			printf(" computing step at iteration %d\n", pmyo->iteration);

	/** next, sort gradient **/
	if (pmyo->verbose) {
		printf("gradient: ");
		for(i = 0; i < n; i++){
			printf("%g ", gradient[i]);
		}
		printf("\n");
	}

	/** gradients_and_indices is the structured array that will only be used as input to the quick sort function 
		Here we prepare the input of qsort by filling in the array
	 **/
	for (i = 0; i < n; i++) {
		gradients_and_indices[i].gradient = gradient[i];
		gradients_and_indices[i].index = i;
	}
	/** Sorting happens here **/
	qsort((void*)gradients_and_indices, n, sizeof(grad_sort_struct), (int(*)(const void*,const void*))compare_grad_components);
	/** Now we output the results in the gradient_sorted array and the sort_index array
		 we use those variables only for convenience and to improve readability
	 **/
	for (i = 0; i < n; i++) {
		gradient_sorted[i] = gradients_and_indices[i].gradient;
		sort_index[i] = gradients_and_indices[i].index;
	}

	/** output some info **/
	if (pmyo->verbose) {
		printf("gradient sorted: ");
		for(i = 0; i < n; i++){
			printf("%g ", gradient_sorted[i]);
		}
		printf("\n");
		printf("indices sorted: ");
		for(i = 0; i < n; i++){
			printf("%d ", sort_index[i]);
		}
		printf("\n");
	}



	/** next, compute direction **/
	num_feasible_y = 0;
	bestcost = 0.0;


	/** Enumerate all possible m like in the pdf **/
	for (m = 0; m < n; m++) {
		/** Initialize y as in (5) **/
		sum = 0.0;
		for (j = 0; j < m; j++) {
			y_sorted[j] = lower[sort_index[j]] - x[sort_index[j]];
			sum += y_sorted[j];
		}
		for (j = m+1; j < n; j++) {
			y_sorted[j] = upper[sort_index[j]] - x[sort_index[j]];
			sum += y_sorted[j];
		}
		y_sorted[m] = -sum; /** use the condition sum of yi = 0 to determine ym**/

		/** check feasibility (lj - xj <= yj <= uj - xj forall j) **/
		y_feasible = 1; /** assume its feasible initially**/
		for (j = 0; j < n; j++) {
			/** if one of the condition is not satisfied, set y_feasible=0 **/
			if (!((lower[sort_index[j]] - x[sort_index[j]] <= y_sorted[j]) &&
					(y_sorted[j] <= upper[sort_index[j]] - x[sort_index[j]]))) {
				y_feasible = 0;
			}
		}
		/** only test if the solution is feasible **/
		if (y_feasible) {
			/** compute the cost = gTy  **/
			cost = 0.0;
			for (j = 0; j < n; j++) {
				cost += gradient_sorted[j] * y_sorted[j];
			}
			/** if its the first time or if the cost is below the best one **/
			if (num_feasible_y == 0 || cost < bestcost) {
				/** found a better solution (or first feasible solution found) **/
				num_feasible_y++;
				/** save the best **/
				bestcost = cost;
				for (j = 0; j < n; j++) {
					y_best_sorted[j] = y_sorted[j];
				}
			}
		}
	}

	if (num_feasible_y == 0) {
		printf("could not find a feasible y, check your bounds.\n");
		retcode = 1;
		goto BACK;
	}
	else {
		if (pmyo->verbose)
			printf("%d feasible(s) directions found\n", num_feasible_y);
		/** Here we print the number of feasible that were found 
		In this case we can determine the number of feasible solutions we expect to see :

		y_j = l_j - x_j if j < m
		y_j = u_j - x_j if j > m

		we have: sum_j y_j = 0 so:
		y_m = - sum_j<m y_j - sum_j>m y_j
		y_m = - sum_j<m l_j + sum_j<m x_j - sum_j>m u_j + sum_j>m x_j
		but sum_j x_j = 1 so:
		y_m + x_m = 1 - sum_j<m l_j - sum_j>m u_j

		In this program the lower bounds are 0 so
		y_m + x_m = 1 - sum_j>m u_j

		y is feasible iff l_j=0 <= y_j + x_j <= u_j for all j
		but this is already the case for all j not equal to m so
		y is feasible iff 0 <= y_m + x_m <= u_m
		y is feasible iff 0 <= 1 - sum_j>m u_j <= u_m
		iff 1 lies in [sum_j>m u_j, sum_j>=m u_j]
		but the partial sums of u_j define an increasing sequence so 1 can lie in at most 1 of the 
		intervals.
		so there can be at most 1 feasible direction in this case.
		 **/



	}

	if (pmyo->verbose)
		printf("best cost: %g\n", bestcost);

	/** set y = y_sorted with the original indices **/
	for (i = 0; i < n; i++) {
		y[sort_index[i]] = y_best_sorted[i];
	}

	/** at this point we have determine the new feasible direction y **/

	if (pmyo->verbose) {
		printf("y* sorted: ");
		for(i = 0; i < n; i++){
			printf("%g ", y_best_sorted[i]);
		}
		printf("\n");


		printf("y*: ");
		for(i = 0; i < n; i++){
			printf("%g ", y[i]);
		}
		printf("\n");
	}


	/** next, compute step size **/

	s = 0.0;

	/** 

	F(x) = lambda x^T Q x - u^T x = lambda x^T V^TFV x + lambda x^T D x - mu^T x
	GradF(x) = 2 lambda V^TFV x + 2 lambda D x - mu
	G(s) = F(x + s y) = s^2 (y^T Q y) + ... with (y^T Q y) > 0
	Since G(s) is a degree two polynomial with positive quadratic coefficient, its a parabola facing upwards
	which implies the existence and unicity of the minimum that we can find by solving G'(s) = 0

	"High level" differentiation of G using a chain rule of matrix calculus:
	(see https://en.wikipedia.org/wiki/Matrix_calculus : Scalar-by-scalar identities With vectors involved)
	G'(s) = (d/ds (x + sy)) . GradF(x + sy) = (d/ds (x + sy))^T * GradF(x + sy)
	G'(s) = y^T * GradF(x + sy)
	G'(s) = y^T GradF(x) + s y^T * (2 lambda V^TFV y + 2 lambda D y)

	We can also show this by expanding F(x + s y) and differentiating each term separately
	but it does not lead immediately to the most compact form.

	we have a linear equation G'(s) = a*s + b = 0 , we compute a and b first then s = -b/a

	we will compute a and b like this :
	a =  2 lambda y^T * (V^TF Vy + D y)
	b = y^T * GradF(x) = y^T * gradient     This is nice because we already computed the gradient

	 **/

	a = 0.0;
	b = 0.0;

	/** compute Vy **/
	for(i = 0; i < f; i++){
		sum = 0;
		for(j = 0; j < n; j++){
			sum += V[i*n + j]*y[j];
		}
		Vy[i] = sum;
	}

	for(j = 0; j < n; j++) {
		/** computes jth entry of VtFVy */
		sum = 0;
		for(i = 0; i < f; i++){
			sum += VtF[j*f + i]*Vy[i];
		}

		/** add (Dy)j == sigma2_j y_j **/
		sum += sigma2[j] * y[j];

		/** multiply by y_j to get the j-th entry of yT * (VTFV y + D y) and add it to a**/
		a += sum * y[j];
	}

	a *= 2.0 * lambda; /** multiply a by its scalor at the end to save some time **/


	/** compute b **/
	sum = 0.0;
	for (j = 0; j < n; j++) {
		sum += y[j] * gradient[j];
	}

	b = sum;

	s = -b/a;

	/** To maintain feasibility we need to constrain s to be in [0, 1]
		Since G(s) is a parabola facing upwards,
		we see that if s is outside [0, 1] we only need to take the closest point of [0, 1] ( 0 or 1 ) to s
		which justifies the following code
	 **/
	if (s > 1.0) {
		s = 1.0;
	}
	if (s < 0.0) {
		s = 0.0;
	}
	/** This second check is however redundant because 
		b = yT * gradient < 0
		a = 2lambda yTDy + yTVTFVy = 2 lambda yT Q y > 0
		so s = -b/a has to be > 0
	 **/

	pmyo->lastcost = pmyo->cost;
	pmyo->cost = bestcost;
	printf("at iter: %d step: %g, cost: %g\n", pmyo->iteration, s, pmyo->cost);

	/** update x **/
	for(j = 0; j < n; j++){
		x[j] += s*y[j];
	}


	BACK:
	return retcode;
}

void myoVtimesy(myo *pmyo, double *y)
{
	/** here y is an n-vector, and we compute V*y and place it in pmyo->Vx **/
	int i, j, f = pmyo->f, n = pmyo->n;
	double sum, *V = pmyo->V;

	for(i = 0; i < f; i++){
		sum = 0;
		for(j = 0; j < n; j++){
			sum += V[i*n + j]*y[j];
		}
		pmyo->Vx[i] = sum;
	}
}

void myo_showx(myo *pmyo, int start, int end)
{
	int j;

	start = start < 0 ? 0 : start;
	end = end >= pmyo->n ? pmyo->n-1 : end;
	/*printf("\n");
	for(j = start; j <= end; j++){
		printf("  x[%d] = %g\n", j, pmyo->x[j]);
	}
	printf("\n");*/

	printf("x: ");
	for(j = start; j <= end; j++){
		printf("%g ", pmyo->x[j]);
	}
	printf("\n");
}

int compare_grad_components(const grad_sort_struct* a, const grad_sort_struct* b) {
	if (a->gradient > b->gradient)
		return -1;
	if (a->gradient < b->gradient)
		return 1;
	return 0;
}

