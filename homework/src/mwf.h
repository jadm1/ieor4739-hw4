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


int MWFMasterThread(int quantity, int numworkers, void **pdatabag,
		int (*masterjobinit)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag),
		int (*masterjobend)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag),
		int (*workerjob)(int jobnumber, pthread_mutex_t *poutputmutex, void* databag, int threadID)
);
int MWFWorkerJobCheckInterrupt(int threadID); /** returns a non zero val if an interrupt signal has been caught (to use in workerjob)**/



#endif
