      function CO(T,n)
	implicit real*8 (A-Z)
	parameter(k=1.380622d-16)!erg K^-1
	if(T.lt.(7.8))then
		CO=0
		return
	endif
	if(dlog10(T).gt.4.5)then
		CO=0
		return
	endif
	T3=T/1000
	Nrat=(3.3d+06)*sqrt(T3)/n
	coeff=4*k*9.7d-8*T*T/2.76/n
	coeff=coeff/(1+Nrat+1.5*sqrt(Nrat))
	CO=coeff
	return
      end
