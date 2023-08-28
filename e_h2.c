/* Feb 21 1997 */
#ifndef _N_BODY
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
	Author: Eraldo Pereira Marinho

		eraldo@iagusp.usp.br

		eraldo@radio.iagusp.usp.br
*/

#define T_rot 85.4 /* kelvin */
#define T_vib 6100 /* kelvin */
#define T_d   52490/* kelvin */
#define N 12
#define J 32
#define para  0
#define ortho 1

#ifndef MINFLOAT
#include <values.h>
#endif

extern double dexp_ (double);

double E_H2 (double T) /* returned as a temperature in kelvins */
{

#if defined (_H_2) || defined (_H_1)

  return 0;

#else

	double tol = (MINFLOAT);
 	register int j;
	register int k=1,e;
  	double E[2], Z[2], z, texp, E_old;
	double control=1;

  if (!T)
  {
	fprintf(stdout,"E_H2: invalid temperature"); fflush(stdout);
	exit (-1);
  }

  for ( j = 0; j < 2; j ++) Z[j]=E[j]=0;

  for ( j = 0; j < J && control>0; j ++)
  {
#undef _VERBOSE
#ifdef _VERBOSE
	k=1-k;
	e=j*(j+1);
	printf("\ne=%i ",e);fflush(stdout);
	texp=1/dexp_(e*T_rot/T);
	printf("texp=%g ",texp);fflush(stdout);
	z=(2*j+1)*texp;
	printf("z=%g ",z);fflush(stdout);
	Z[k]+=z;
	printf("Z[%i]=%g ",k,Z[k]);fflush(stdout);
	E_old=E[k];
	printf("E[%i]=%g ",k,E[k]);fflush(stdout);
	E[k]+=z*e;
	printf("E[%i]=%g ",k,E[k]);fflush(stdout);
	if (j>4)
          control=
	    fabs(E[k]-E_old)-tol*(E[k]+E_old);
	    printf("control=%g\n",control);fflush(stdout);
#else
	k	=	1-k;
	e	=	j*(j+1);
	texp	=	1/dexp_(e*T_rot/T);
	z	=	(2*j+1)*texp;
	Z[k]	+=	z;
	E_old   =	E[k];
	E[k]	+=	z*e;
	if (j>4)
          control =
	    fabs(E[k]-E_old)-tol*(E[k]+E_old);
#endif
  }

  return
    T_rot*(E[para]+3*E[ortho])/(Z[para]+3*Z[ortho])
    +
    T_vib/(dexp_(T_vib/T)-1);

#endif

}

#endif
