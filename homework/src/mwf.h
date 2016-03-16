#ifndef MWF_H
#define MWF_H


#define WAITING 100
#define WORKING 101
#define PREANYTHING 102
#define DONEWITHWORK 103


#define QUIT 200
#define WORK 201
#define STANDBY 202
#define INTERRUPT 203

/**
 * Master worker framework header file
 * **/


typedef struct WorkerBag {
	int ID; /** worker thread ID **/
	int status; /** status code **/
	int command; /** command code **/
	int jobnumber;
	int itercount;
	pthread_mutex_t *psynchro; /** mutex pointer for communication with the master thread **/
	pthread_mutex_t *poutputmutex; /** mutex pointer for outputing text to the console**/
	int (*masterjobinit)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag);
	int (*masterjobend)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag);
	int (*workerjob)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag, int threadID);
	void *databag;
} WorkerBag;


int MWFMasterThread(int quantity, int numworkers, void **pdatabag,
		int (*masterjobinit)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag),
		int (*masterjobend)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag),
		int (*workerjob)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag, int threadID)
);
int MWFWorkerJobCheckInterrupt();



#endif
