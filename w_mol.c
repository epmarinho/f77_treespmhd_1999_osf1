#ifndef _N_BODY
#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>

#  ifdef _LINUX_
float w_mol__
#  else
float w_mol_
#  endif
	(float *X,float *Y,float *y)
{
  return .5*(*X*(*y+1)+*Y*.5);
}
#endif
