
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "utilities.h"

/**
 * Safe freeing
 * free an address and set it to NULL to prevent double freeing
 */
void UTLFree(void **paddress)
{
	void *address = *paddress;

	if (address == NULL) goto BACK;

	/**printf("freeing array at %p\n", address);**/
	free(address);
	address = NULL; /** prevents double freeing **/

	BACK:
	*paddress = address;
}


/** print vector **/
void UTLShowVector(int n, double *vector)
{
	int j;

	for (j = 0; j < n; j++){
		printf(" %g", vector[j]);
	}
	printf("\n");
}


char *Ggettimestamp(void)
{
  time_t timestamp;

  timestamp = time(0);
  return (char *) ctime(&timestamp);
}




#ifdef WIN32
/* Reentrant random function frm POSIX.1c.
   Copyright (C) 1996, 1999 Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Ulrich Drepper <drepper@cygnus.com <mailto:drepper@cygnus.com>>, 1996.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */
   
/* This algorithm is mentioned in the ISO C standard, here extended
   for 32 bits.  */
int rand_r (unsigned int *seed)
{
    unsigned int next = *seed;
    int result;

    next *= 1103515245;
    next += 12345;
    result = (unsigned int) (next / 65536) % 2048;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next / 65536) % 1024;

    next *= 1103515245;
    next += 12345;
    result <<= 10;
    result ^= (unsigned int) (next / 65536) % 1024;

    *seed = next;

    return result;
}

#endif

