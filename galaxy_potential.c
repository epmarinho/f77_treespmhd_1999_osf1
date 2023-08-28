#include <math.h>
#include <stdlib.h>
#include "special.h"

#define CONV_PHI .0006536842725625
#define CONV_G   .006536842725625


double bessi0(double);
double bessk0(double);
double bessi1(double);
double bessk1(double);

/* recomended Rcore < 0.1 kpc */
float phi_galax_(float x[3], float g[3], float *R)
{
   double Rcore=*R;
   double eb=1.0/4.81;
   double mb=53500.00;
   double md=16000.00;
   double md2=6000.00;
   double ed=1.0/2.60;
   double ed2=1.0/4.50;
   double sigma0=300.00;
   double p=0.90;
   double a1=12250.0;
   double a2=-1475.0;
   double b1=4.20;
   double b2=2.25;
   double delt=0.5;
   double rp,anel,vgas,mi;
   double vbojo,ad,vdisk,vdisk2,rot,ad2;
   double phig,phib,phid,phidz,phit,r,r3;
   double pi=4*atan(1.0);
   double g_aux;
   double pvbojo;
   
   r=x[0]*x[0] + x[1]*x[1] + Rcore;
   rp = r * (p*p)*0.5;
   r3= r + x[2]*x[2];
   r = sqrt(r); r3 = sqrt(r3);

   mi=647*exp(-r*ed)+203*exp(-r*ed2)+11.00;
   phig=-(pow(pi,1.5))*4.3*sigma0*(1.0/p)*exp(-rp)*bessi0(rp);
   phig=-(4/3)*pi*4.3*a1*(((1.0)/(b1*pow((b1*b1+r*r),0.5)))-
	 ((0.5*b1)/(pow((b1*b1+r*r),1.5))))+phig;
   phig=-(4/3)*pi*4.3*a2*(((1.0)/(b2*pow((b2*b2+r*r),0.5)))-
	 ((0.5*b2*r*r)/(pow((b2*b2+r*r),1.5))))+phig;
 
  
   /* componente bojo */
   phib=-mb/(r3+eb);

   
   /* componente disco - md = 2.pi.G.mio*/
   ad=ed*r/2.0;
   phid=-((md/2.0)*r*(
	 (bessi0(ad)*bessk1(ad))-(bessi1(ad)*bessk0(ad))
	 ));
   ad2=ed2*r/2;
   phid=-((md2/2.0)*r*(
	  (bessi0(ad2)*bessk1(ad2))-(bessi1(ad2)*bessk0(ad2))
	  ))+phid;
 

  /* componente z do potencial do disco */
   phidz=2*pi*4.3*mi*delt*log(2*cosh(x[2]/delt));
  

  /* potencial total */
   phit=phib+phid+phig+phidz;


   /* calculo da velocidade devido ao gas*/
   rp=(pow((r*p),2.0))/2.0;
   anel=(pow(pi,1.5))*sigma0*r*r*p*exp(-rp)*((bessi0(rp)
	 -bessi1(rp)))*4.3;
   anel=(4/3)*pi*a1*(((r*r)/(b1*pow((b1*b1+r*r),1.5)))-
	 ((1.5*b1*r*r)/(pow((b1*b1+r*r),2.5))))*4.3+anel;
   anel=(4/3)*pi*a2*(((r*r)/(b2*pow((b2*b2+r*r),1.5)))-
	 ((1.5*b2*r*r)/(pow((b2*b2+r*r),2.5))))*4.3+anel;
   vgas=pow((fabs(anel)),0.5);

   /* componente bojo */
   vbojo=pow(((mb*r3)/(pow((r3+eb),2.0))),0.5);

   /* componente disco */
   vdisk=pow((md*r*ad*((bessi0(ad)*bessk0(ad))-
	 (bessi1(ad)*bessk1(ad)))),0.5);
   vdisk2=pow((md2*r*ad2*((bessi0(ad2)*bessk0(ad2))-
	  (bessi1(ad2)*bessk1(ad2)))),0.5);

   /* velocidade total */
   /* rot=pow(((pow(vdisk,2.0))+ (pow(vdisk2,2.0))+ (pow(vgas,2.0))+
       (pow(vbojo,2.0))),0.5); */ 

   pvbojo = pow(vbojo,2.0);

   g_aux = CONV_G * (1/(r*r))*(pow(vdisk,2.0)+pow(vdisk2,2.0)+pow(vgas,2.0))+
	 (x[0]/(r3*r3))*pvbojo;

   g[0]= - x[0] * g_aux;

   g[1]= - x[1] * g_aux;

   g[2]= - CONV_G * ((x[2]/(r3*r3))*pvbojo + 2*pi*4.3*mi*tanh(x[2]/delt));

  return CONV_PHI * phit;
}
