#include <stdio.h>
#include <stdlib.h>
#include "myopt.h"

int MYOcreatemyo(myo **pmyo);

int main(int argc, char *argv[])
{
  int code = 0;
  myo *pmyo = NULL;

  if(argc < 2){ 
    printf(" usage: myo filename\n");
    code = 1; goto BACK;
  }

  printf("%s\n", Ggettimestamp());

  if( (code = myoGetmyoFromFile(&pmyo, argv[1])) ) goto BACK;

  code = myoalgo(pmyo);

  myokillmyo(&pmyo);
	 
BACK:
  return code;
}
