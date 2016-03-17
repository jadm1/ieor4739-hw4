/**
 * Master worker framework source file
 * **/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>
#include "utilities.h"
#include "mwf.h"

/**
 * Local worker bag structure
 */
typedef struct WorkerBag {
	int ID; /** worker thread ID **/
	int status; /** status code **/
	int command; /** command code **/
	int jobnumber;
	pthread_mutex_t *psynchro; /** mutex pointer for communication with the master thread **/
	pthread_mutex_t *poutputmutex; /** mutex pointer for outputing text to the console**/
	int (*masterjobinit)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag);
	int (*masterjobend)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag);
	int (*workerjob)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag, int threadID);
	void *databag;
} WorkerBag;


static WorkerBag **pwbagproxy = NULL;
static int numworkersproxy = 0;
static char *deadstatus = NULL;
static int activeworkers = 0;

#ifdef WIN32
#define sigset signal
#endif
void (*sigset(int sig, void (*disp)(int)))(int);
void handlesigint(int i);

void* WorkerThread(void *voidedwbag);



void handlesigint(int signal)
{
	int j;
	printf("yo, what's happening\n");
	for(j = 0; j < numworkersproxy; j++){
		deadstatus[j] = 1;
	}
	/** brutal **/
}

int MWFWorkerJobCheckInterrupt(int threadID) {
	int interrupting = 0;
	WorkerBag **pwbag = pwbagproxy;
	WorkerBag *wbag;

	if (pwbag == NULL)
		return 0;

	wbag = pwbag[threadID];

	pthread_mutex_lock(wbag->psynchro);

	interrupting = 0;
	if (wbag->command == INTERRUPT || wbag->command == QUIT){
		interrupting = 1;
	}

	pthread_mutex_unlock(wbag->psynchro);

	return interrupting;
}


int MWFMasterThread(int quantity, int numworkers, void **pdatabag,
		int (*masterjobinit)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag),
		int (*masterjobend)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag),
		int (*workerjob)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag, int threadID)
) {
	int retcode = 0;
	int j;
	WorkerBag **pwbag = NULL;
	WorkerBag *wbag = NULL;
	int initialruns, scheduledjobs;
	int theworker;
	char gotone;
	pthread_t *pthread;
	pthread_mutex_t outputmutex;
	pthread_mutex_t *psyncmutex;

	/**
	 * Master Thread
	 */

	if ( numworkers > quantity ){
		numworkers = quantity;
		printf(" --> reset workers to %d\n", numworkers);
	}

	deadstatus = (char *) calloc(numworkers, sizeof(char));
	sigset(SIGINT, &handlesigint);

	pthread_mutex_init(&outputmutex, NULL); /** common to everybody **/

	psyncmutex = (pthread_mutex_t *)calloc(numworkers, sizeof(pthread_mutex_t));
	if(!psyncmutex){
		printf("could not create mutex array\n"); retcode = NOMEMORY; goto BACK;
	}

	for(j = 0; j < numworkers; j++)
		pthread_mutex_init(&psyncmutex[j], NULL);

	pwbag = (WorkerBag **)calloc(numworkers, sizeof(WorkerBag *));
	if(!pwbag){
		printf("could not create mwf bag array\n"); retcode = NOMEMORY; goto BACK;
	}
	for (j = 0; j < numworkers; j++) {
		pwbag[j] = (WorkerBag *)calloc(1, sizeof(WorkerBag));
	}

	pwbagproxy = pwbag;
	numworkersproxy = numworkers;

	pthread = (pthread_t *)calloc(numworkers, sizeof(pthread_t));
	if (pthread == NULL) {
		printf("could not create thread array\n"); retcode = NOMEMORY; goto BACK;
	}

	/**
	 * Initialize worker bags
	 */
	for (j = 0; j < numworkers; j++) {
		wbag = pwbag[j];
		wbag->ID = j;
		wbag->command = STANDBY;
		wbag->status = PREANYTHING;
		wbag->psynchro = &psyncmutex[j];
		wbag->poutputmutex = &outputmutex;
		wbag->databag = pdatabag[j];
		wbag->masterjobinit = masterjobinit;
		wbag->masterjobend = masterjobend;
		wbag->workerjob = workerjob;
	}

	/**
	 * Launch worker threads
	 */

	for(j = 0; j < numworkers; j++) {
		printf("Launching thread for worker %d\n", j);
		pthread_create(&pthread[j], NULL, &WorkerThread, (void *) pwbag[j]);
	}

	initialruns = numworkers;
	if (initialruns > quantity)
		initialruns = quantity;

	for(theworker = 0; theworker < initialruns; theworker++) {
		wbag = pwbag[theworker];

		/**
		 * Master job init here
		 */
		if (masterjobinit != NULL)
			masterjobinit(wbag->jobnumber, &outputmutex, (void*)wbag->databag);

		pthread_mutex_lock(&outputmutex);
		printf("*****master:  worker %d will run job %d\n", theworker, theworker);
		pthread_mutex_unlock(&outputmutex);

		/** tell the worker to work **/
		pthread_mutex_lock(&psyncmutex[theworker]);
		wbag->command = WORK;
		wbag->status = WORKING;
		wbag->jobnumber = theworker;
		pthread_mutex_unlock(&psyncmutex[theworker]);

	}
	scheduledjobs = activeworkers = initialruns;

	while (activeworkers > 0) {
		/** check the workers' status **/
		gotone = 0;
		for(theworker = 0; theworker < numworkers; theworker++){

			pthread_mutex_lock(&psyncmutex[theworker]);
			wbag = pwbag[theworker];
			if(wbag->status == DONEWITHWORK){

				pthread_mutex_lock(&outputmutex);
				printf("master:  worker %d is done with job %d\n", wbag->ID, wbag->jobnumber);
				pthread_mutex_unlock(&outputmutex);

				/**
				 * Master job end here
				 */
				if (masterjobend != NULL)
					masterjobend(wbag->jobnumber, &outputmutex, (void*)wbag->databag);

				if(scheduledjobs >= quantity){
					/** tell worker to quit **/
					/**pthread_mutex_lock(&outputmutex);
						printf("master: telling worker %d to quit\n", wbag->ID);
						pthread_mutex_unlock(&outputmutex);**/
					wbag->command = QUIT;
					wbag->status = QUIT;
					--activeworkers;
				}
				else {
					gotone = 1;
				}
			}
			else if(wbag->status == PREANYTHING) {
				/**pthread_mutex_lock(&outputmutex);
					printf("master:  worker %d is available\n", theworker);
					pthread_mutex_unlock(&outputmutex);**/
				gotone = 1;
			}
			/**else if( (wbag->status == WORKING) && (wbag->itercount > 100000)){
				wbag->command = INTERRUPT;
				pthread_mutex_lock(&outputmutex);
				printf("master: telling worker %d to interrupt\n", wbag->ID);
				pthread_mutex_unlock(&outputmutex);
			} not keeping track of iterations anymore**/
			else if(deadstatus[theworker]){
				/**pthread_mutex_lock(&outputmutex);
					printf("master: telling worker %d to quit\n", wbag->ID);
					pthread_mutex_unlock(&outputmutex);**/
				wbag->command = QUIT;
				wbag->status = QUIT;
				--activeworkers;
				/** and let's make sure we don't do it again **/
				deadstatus[theworker] = 0;
			}
			pthread_mutex_unlock(&psyncmutex[theworker]);
			if(gotone) break;
			usleep(10000);

		}
		/** at this point we have run through all workers **/

		if(gotone){
			/** if we are here, "theworker" can work **/
			wbag = pwbag[theworker];

			if (masterjobinit != NULL)
				masterjobinit(wbag->jobnumber, &outputmutex, (void*)wbag->databag);

			pthread_mutex_lock(&outputmutex);
			printf("master:  worker %d will run job %d\n", theworker, scheduledjobs);
			pthread_mutex_unlock(&outputmutex);

			/** tell the worker to work **/
			pthread_mutex_lock(&psyncmutex[theworker]);
			wbag->command = WORK;
			wbag->status = WORKING;
			wbag->jobnumber = scheduledjobs;
			pthread_mutex_unlock(&psyncmutex[theworker]);

			++scheduledjobs;
		}
	}



	/*  pthread_mutex_lock(&psynchro_array[theworker]);
	  pbag->command = QUIT;
	  pthread_mutex_unlock(&psynchro_array[theworker]);*/

	pthread_mutex_lock(&outputmutex);
	printf("master:  done with loop\n");
	pthread_mutex_unlock(&outputmutex);


	/**
	 * Waiting for worker threads to terminate
	 */
	for(j = 0; j < numworkers; j++){
		pthread_join(pthread[j], NULL);
		pthread_mutex_lock(&outputmutex);
		printf("master: joined with thread %d\n", j);
		pthread_mutex_unlock(&outputmutex);
	}

	/**
	 * Free worker bags
	 */
	for(j = 0; j < numworkers; j++){
		wbag = pwbag[j];
		UTLFree((void**)&pwbag[j]);
	}
	UTLFree((void**)&pwbag);



	BACK:
	return retcode;
}




void* WorkerThread(void *voidedwbag)
{
	WorkerBag *wbag = (WorkerBag *)voidedwbag;
	int letsgo = 0, forcedquit = 0;
	int waitcount;
	int (*workerjob)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag, int threadID);
	workerjob = wbag->workerjob;

	pthread_mutex_lock(wbag->poutputmutex);
	printf("ID %d starts\n", wbag->ID);
	pthread_mutex_unlock(wbag->poutputmutex);




	for(;;){
		/**pthread_mutex_lock(wbag->poutputmutex);
		printf(" ID %d in big loop\n", wbag->ID);
		pthread_mutex_unlock(wbag->poutputmutex);**/

		letsgo = 0;
		waitcount = 0;
		while(letsgo == 0){
			/** wait until WORK signal **/
			usleep(10000);

			pthread_mutex_lock(wbag->psynchro);
			if(wbag->command == WORK){
				letsgo = 1;
			}
			else if(wbag->command == QUIT)
				letsgo = 2;
			pthread_mutex_unlock(wbag->psynchro);

			if (letsgo == 2)
				goto BACK;

			if(0 == waitcount%20){
				/**pthread_mutex_lock(wbag->poutputmutex);
				printf("ID %d bag %p: wait %d for signal; right now have %d\n", wbag->ID, (void *) wbag, waitcount, wbag->command);
				pthread_mutex_unlock(wbag->poutputmutex);**/
			}
			++waitcount;

		}

		/**pthread_mutex_lock(wbag->poutputmutex);
		printf("ID %d: got signal to start working\n", wbag->ID);
		pthread_mutex_unlock(wbag->poutputmutex);**/


		/**
		 * worker job here
		 */

		if (workerjob != NULL)
			workerjob(wbag->jobnumber, wbag->poutputmutex, (void*)wbag->databag, wbag->ID);




		/** first, let's check if we have been told to quit **/
		pthread_mutex_lock(wbag->psynchro);
		if(wbag->command == QUIT)
			forcedquit = 1;
		pthread_mutex_unlock(wbag->psynchro);

		if(forcedquit)
			break;

		pthread_mutex_lock(wbag->psynchro);
		wbag->status = DONEWITHWORK;
		wbag->command = STANDBY;
		pthread_mutex_unlock(wbag->psynchro);
	}

	BACK:
	pthread_mutex_lock(wbag->poutputmutex);
	printf(" ID %d quitting\n", wbag->ID);
	pthread_mutex_unlock(wbag->poutputmutex);

	return (void *)&wbag->ID;
}




