/*
 * utilities.h
 *
 *  Created on: Mar 16, 2016
 *      Author: root
 */

#ifndef UTILITIES_H
#define UTILITIES_H

#define NOMEMORY 100

void UTLFree(void **paddress);
void UTLShowVector(int n, double *vector);
char *Ggettimestamp(void);
#ifdef WIN32
int rand_r (unsigned int *seed);
#endif


#endif
