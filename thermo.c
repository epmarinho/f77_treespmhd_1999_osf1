#ifndef _N_BODY

#include <stdlib.h>
#include <math.h>
#undef _VERBOSE
#ifdef _VERBOSE
#include <stdio.h>
#endif

#ifndef MAXDOUBLE
#include <values.h>
#endif

#define	D0	5.1955e+4	/* Kelvin */
#define	Ro	8.154319e-3	/* [e/m] K^-1 */
#define X       0.75
#define Y       0.25

#if ! ( defined(_TUBE) || defined(_H_1) || defined(_H_2) )

double dexp_(double x)
{
    double S=0;
    int n,N=24;
    int sgn=x<0?1:0;
    x=fabs(x);
    for(n=N; n>0 && x*S/(n+1)<MAXDOUBLE; n--) S=1+x*(S/n);
    if(sgn)S=1/S;
	return S;
}

float fexp_(float x)
{
	 double S=0;
	 int n,N=24;
	 int sgn=x<0?1:0;
	 x=fabs(x);
	 for(n=N; n>0 && x*S/(n+1)<MAXDOUBLE; n--) S=1+x*(S/n);
	 if(sgn)S=1/S;
	return (float)S;
}

/*
	Author: Eraldo Pereira Marinho

		eraldo@iagusp.usp.br

		eraldo@radio.iagusp.usp.br
*/


#ifndef MINFLOAT
#include <values.h>
#endif

#define u_therm	0.956132e+10
#define	u_dens	1.504929567e-20 /* g cm^-3 */

extern double dexp_(double);

extern double E_H2(double);

#ifdef _LINUX_
void y_weight__
#else
void y_weight_
#endif
(float *rho, float *T, float *y)
{
#ifdef _H_2
    *y = 0;
#elif _H_1
    *y = 1;
#else
    double y_ = 0;
    void y_weight( double, double, double*);
    y_weight( *rho, *T, &y_);
    *y = y_;
#endif
}

void y_weight(rho, T, y) /* Warning: CGS units! */
double rho, T, *y;
{  double A, e;
  e = D0/T;
  A = dexp_(e)*rho;
  A = 2.11/A;
  *y=0;
  if (A>1)
  {
    if (A<.03125e9) *y = .5*A * ( sqrt(1+4/A) - 1 );
    else *y = 1 - 2/A * ( 1-6/A );
  }
  else *y = sqrt(A*(.25*A+1)) - .5*A;
}

#endif
double	rhoCGS;

float temperature_(
	float *u,
	float *rho,
	float *P,
	float *over_mu,
	float *yy,
	float *freedom
)
{
	double	T;

#if ! ( defined(_TUBE) || defined(_H_1) || defined(_H_2) )
  double	control = 1.;
  double	tol_T   = (MINFLOAT);
  double	old_T;
  double	EH2;
  register      itc     = 0;
#endif

#ifdef _H_2
#  define	w_m	0.5
#  define	y	0.0
  *freedom = 5;
#elif  _H_1
#  define	w_m	1.0
#  define	y	1.0
  *freedom = 3;
#else
  double	y = *yy>0?*yy:0; /* November 8th 1996 */
  double	w_m=.5*(X*(1+y)+.5*Y);
  *freedom = 4;
#endif

  /* guess an initial temperature */
  T=(*u)/(.5*(*freedom)*w_m*Ro);
#if !( defined (_H_2) || defined (_H_1) )

  rhoCGS = *rho*u_dens;

  while (control>0) {
#ifdef _VERBOSE
    old_T=T;printf("T=%g ",T);fflush(stdout);
    y_weight(rhoCGS*X,T,&y);printf("y=%g ",y);fflush(stdout);
    w_m=.5*(X*(1+y)+.5*Y);printf("m=%g ",1/w_m);fflush(stdout);
    EH2=E_H2(T);printf("E=%g ",EH2);fflush(stdout);
    T=(*u/Ro-.5*X*(y*D0+(1-y)*EH2))/(w_m*1.5);printf("T=%g ",T);fflush(stdout);
#else
    old_T = T;
    y_weight(rhoCGS*X,T,&y);
    w_m = .5*(X*(1+y)+.5*Y);
    EH2 = E_H2(T);
    T	= (*u/Ro-.5*X*(y*D0+(1-y)*EH2))/(w_m*1.5);
#endif
    if (T<0) T = .5 * old_T;
    control = fabs(T-old_T)-tol_T*(T+old_T);
    itc++;
    if (itc>0x80) {
      T=.5*(T+old_T);
      if (itc>0x100) control = -1;
    }
  }

  /*
	*P is returned as P/rho^2:

	P = rho * R*T/mu

	P/rho^2 = R*T/(mu*rho)

  */

  *P = (float) ( (Ro*T*w_m)/(*rho) );
  *freedom = (float) ( (*u)/(.5*T*w_m*Ro) );

#else

  *P = (float) ( (2/(*freedom))*(*u)/(*rho) );

#endif

  *yy = (float) y;
  *over_mu = (float) w_m;
  return (float) T;

}

#endif
