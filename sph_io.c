#include <stdio.h>
#include <stdlib.h>

void prl_(int *itab,char *strng)
{
	int i;
	fflush(stdout);
	for(i=0;i<*itab;i++)printf(" ");
	printf("%s",strng);
	fflush(stdout);
}
