      subroutine cooling_H_He(T,n,LAMBDA)
	implicit real*8 (A-Z)
	parameter(T1=5.23,T2=6.2,z=0)
!
!	[rho] ~= 1.51e-19 g cm^-3
!	[n]   ~= 9.09e+04 cm^-3
!
	logT = log10(T)
	if(logT.lt.(3.8))then
	  LAMBDA=0
	  return
	endif
	texp1 = (T1-logT)
	texp1 = texp1*texp1
	texp1 = -.1-1.88*(texp1*texp1)
	result= 0
	call exponentialize(texp1,result)
	texp2 = 0
	if(logT.le.T2)then
	    texp2 = (T2-logT)
	    texp2 = texp2*texp2
	    texp2 = 0.2*(texp2*texp2)
	endif
	texp2 = -1.7-texp2
	call exponentialize(texp2,result)
	result = n*n * result * 1.0d-21
	LAMBDA = result
	! adding up the inverse Comptom cooling:
	ne=n
	result=0
	if(logT.gt.4)result=5.41d-36*ne*T
	! for primordial times:
	!if(logT.gt.4)result=5.41d-36(1+z)**4*ne*T
	LAMBDA = result+LAMBDA 
	return
      end
!
!
      subroutine exponentialize(xexp,result)
	implicit real*8 (A-Z)
	if(xexp.gt.0)then
	  if(xexp.lt.38)then
	    result = result + 10.0**xexp
	  else
	    result = result + 10.0**38
	  endif
	else
	  if(xexp.gt.-38)result = result + 10.0**xexp
	endif
	return
      end
!
!
      subroutine cooling_rate_2(n,T,LAMBDA)
	real LAMBDA
!
	if(T.lt.300)then
		result=0
	else
	if(T.lt.2000)then
		result=2.238e-32*(T*T)
	else
	if(T.lt.8000)then
		result=1.0012e-30*(T**1.5)
	else
	if(T.lt.1e+5)then
		result=4.6240e-36*(T**2.867)
	else
	if(T.lt.4e+7)then
		result=1.78e-18*(T**-.65)
	else
		result=3.2217e-27*(T**.5)
	endif
	endif
	endif
	endif
	endif
	LAMBDA=result
	return
      end
