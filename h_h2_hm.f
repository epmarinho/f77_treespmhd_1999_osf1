      function g_10_H(T)
	implicit real*8 (a-z)
	parameter(c=1.0d-12)
	T13=1000/T
	if(T13.gt.300)then
	  g_10_H=0
	  return
	endif
	g_10_H=c*dsqrt(T)/dexp(T13)
	return
      end

      function g_20_H(T)
	implicit real*8 (a-z)
	parameter(c=1.6d-12)
	aux=400/T
	aux=aux*aux
	if(aux.gt.300)then
	  g_20_H=0
	  return
	endif
	g_20_H=c*dsqrt(T)/dexp(aux)
	return
      end

      function g_10_H2(T)
	implicit real*8 (a-z)
	parameter(c=1.4d-12)
	parameter(d=1.2d3,f=1.2d4)
	e=f/(T+d)
	g_10_H2=c*sqrt(T)/dexp(e)
	return
      end

      function g_20_H2(T)
	implicit real*8 (a-z)
	g_20_H2=0.0d+00
	return
      end

      function g_H(J,T)
	implicit real*8 (a-z)
	integer J
	parameter(a=1.0d-11,b=6.0d+01,c=1.0d-12)
	parameter(d=.33,e=0.9,f=-3.5)
	T3=T/1000
	t34=1/T3
	t34=t34*t34
	t34=t34*t34
	g=(J+f)/e
	g=g*g
	coeff=a*dsqrt(T3)/(1+b*t34)+c*T3
	coeff=coeff*( d + e/dexp(g) ) !cm^3/sec
	g_H=coeff
	return
      end

      function g_H2(J,T)
	implicit real*8 (a-z)
	integer J
	T3=T/(1d+3)
	coeff=(3.3d-12)+(6.6d-12)*T3
	coeff=coeff*(.276*J*J/dexp((J/3.18d0)**1.7))!cm^3/sec
	g_H2=coeff! cm^3/sec
	return
      end

      function L_r(T)
	implicit real*8 (a-z)
	T3=T/1d+3
	T13=.13/T3
	T13=T13*T13*T13
	if(T13.gt.300)then
	  L_r=0
	  return
	endif
	coeff=(9.5d-22*(T3**3.76))/(1+.12*T3**2.1)
	coeff=coeff/dexp((.13/T3)**3)+3d-24/dexp(.51/T3)
	L_r=coeff! erg/sec
	return
      end

      function L_v(T)
	implicit real*8 (a-z)
	T3=T/1d+3
	T13=5.86/T3
	if(T13.gt.300)then
	  L_v=0
	  return
	endif
	coeff=6.7d-19/dexp(5.86/T3)
	coeff=coeff+1.6d-18/dexp(11.7/T3)
	L_v=coeff! erg/sec
	return
      end

      function L_low_r_H(T)
	implicit real*8 (a-z)
	parameter(E0=0,E1=.8583087956d+02)
	E2=3*E1
	E3=6*E1
	g_2 = g_H(2,T)
	g_3 = g_H(3,T)
        
	coeff=0
	if(((E2-E0)/T).lt.300)then
	  coeff=.25*(5*g_2*(E2-E0)/dexp((E2-E0)/T))
	endif
	if(((E3-E1)/T).lt.300)then
	  coeff=coeff+.75*(7./3*g_3*(E3-E1)/dexp((E3-E1)/T))
	endif
	L_low_r_H=coeff! erg cm^3/sec
	return
      end

      function L_low_v_H(T)
	implicit real*8 (a-z)
	parameter(E10=5.86d+3)
	E20=2*E10
	g_10 = g_10_H(T)
	g_20 = g_20_H(T)
	coeff=0
	if((E10/T).lt.300)then
	  coeff=g_10*E10/dexp(E10/T)
	endif
	if((E20/T).lt.300)then
	  coeff=coeff+g_20*E20/dexp(E20/T)
	endif
	L_low_v_H=coeff! erg cm^3/sec
	return
      end

      function L_low_r_H2(T)
	implicit real*8 (a-z)
	parameter(E0=0,E1=.8583087956d+02)
	E2=3*E1
	E3=6*E1
	g_2 = g_H2(2,T)
	g_3 = g_H2(3,T)
	coeff=0
	if(((E2-E0)/T).lt.300)then
	  coeff=2.5d-1*(5*g_2*(E2-E0)/dexp((E2-E0)/T))
	endif
	if(((E3-E1)/T).lt.300)then
	  coeff=coeff+7.5d-1*(7./3*g_3*(E3-E1)/dexp((E3-E1)/T))
	endif
	L_low_r_H2=coeff! erg cm^3/sec
	return
      end

      function L_low_v_H2(T)
	implicit real*8 (a-z)
	parameter(E10=5860)
	E20=2*E10
	g_10 = g_10_H2(T)
	g_20 = g_20_H2(T)
	coeff=0
	if((E10/T).lt.300)then
	  coeff=g_10*E10/dexp(E10/T)
	endif
	if((E20/T).lt.300)then
	  coeff=coeff+g_20*E20/dexp(E20/T)
	endif
	L_low_v_H2=coeff! erg cm^3/sec
	return
      end

      subroutine cooling_H_H2(T,n_H,n_H2,L_H,L_H2)
	implicit real*8 (a-z)

	lT=dlog10(T)
	if(lT.gt.4.5)then
	  L_H=0
	  L_H2=0
	  return
	endif

	Lrot=L_r(T)
	Lvib=L_v(T)
	Llow_r=L_low_r_H(T)
	Llow_v=L_low_v_H(T)

	! rotational due to H collisions:
	L_H_r=Lrot*Llow_r
	if(L_H_r.gt.0)then
	    if(Lrot.gt.0)L_H_r=L_H_r/(n_H*Llow_r+Lrot)!ergs cm^3/sec
	endif
	
	! vibrational due to H collisions:
	L_H_v=Lvib*Llow_v
	if(L_H_v.gt.0)then
	    if(Lvib.gt.0)L_H_v=L_H_v/(n_H*Llow_v+Lvib)!ergs cm^3/sec
	endif

	L_H=L_H_r+L_H_v

	Llow_r=L_low_r_H2(T)
	Llow_v=L_low_v_H2(T)

	! rotational due to H2 collisions:
	L_H2_r=Lrot*Llow_r
	if(L_H2_r.gt.0)then
	    if(Lrot.gt.0)L_H2_r=L_H2_r/(n_H2*Llow_r+Lrot)!ergs cm^3/sec
	endif
	
	! vibrational due to H2 collisions:
	L_H2_v=Lvib*Llow_v
	if(L_H2_v.gt.0)then
	    if(Lvib.gt.0)L_H2_v=L_H2_v/(n_H2*Llow_v+Lvib)!ergs cm^3/sec
	endif

	L_H2=L_H2_r+L_H2_v

	return
      end
