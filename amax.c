#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
float amax_(float *x, float *y)
{
  return (*x>*y ? *x : *y);
}

int max_(int *x, int *y)
{
  return (*x>*y ? *x : *y);
}
