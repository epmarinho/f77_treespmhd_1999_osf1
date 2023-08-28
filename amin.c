#include <math.h>
#include <stdlib.h>

float amin_(float *x, float *y)
{
  return (*x<*y ? *x : *y);
}

int min_(int *x, int *y)
{
  return (*x<*y ? *x : *y);
}
