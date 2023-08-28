      function cooling(n,T)
	implicit real (a-z)
	common/H_H2_CO/ Lamb(0:1023,0:63),Hy(0:1023,0:63)
	integer i_n,i_T
	if(n.eq.0)then
	  cooling=0
	  return
	endif
	if(T.eq.0)then
	  cooling=0
	  return
	endif
	lgn=log10(n)
	lgT=log10(T)
	!print*,lgn,lgT
	x_n=(lgn+2)*64/13
	x_T=lgT*1024/7
	!print*,x_n,x_T
	i_n=int(x_n)
	if(i_n.lt.0)then
	  cooling=0
	  return
	endif
	i_T=int(x_T)
	if(i_T.lt.0)then
	  cooling=0
	  return
	endif
	if(i_n.gt.63)i_n=63
	if(i_T.gt.1023)i_T=1023
	!print*,i_n,i_T
	L0=Lamb(i_T,i_n)
	i_n1=i_n
	if(i_n.lt.63)i_n1=i_n+1
	D_L_n=Lamb(i_T,i_n1)-L0
	i_T1=i_T
	if(i_T.lt.1023)i_T1=i_T+1
	D_L_T=Lamb(i_T1,i_n)-L0
	D_n=x_n-i_n
	D_T=x_T-i_T
	D_L=D_L_n*D_n+D_L_T*D_T
	L=L0+D_L
	L=10**L
	!print*,L
	!warning: positively returned value.
	! set efficience to be 25%
	!*** L=.25*L
	cooling=L
	return
      end
!
!
      subroutine Assembl_Table
	implicit real (a-z)
	common/H_H2_CO/ Lamb(0:1023,0:63),Hy(0:1023,0:63)
	integer i_n,i_T
	open(1,file='H+H2+CO.dat',status='old')
	print*,'Reading the cooling table.'
	print*,'Please, wait.'
	do i_n=0,63
	  do i_T=0,1023
	    read(1,*)x,x,Lamb(i_T,i_n),Hy(i_T,i_n)
	  end do
	end do
	close(1)
	print*,'Ready.'
	print*,'Writing the binary cooling table.'
	print*,'Please, wait.'
	open(1,file='H+H2+CO.bin'
     &	,status='unknown',form='unformatted')
	write(1)Lamb
	close(1)
	print*,'Ready.'
	return
      end
!
      function H_mol(n,T)
	implicit real (a-z)
	common/H_H2_CO/ Lamb(0:1023,0:63),Hy(0:1023,0:63)
	integer i_n,i_T
	if(n.eq.0)then
	  H_mol=0
	  return
	endif
	if(T.eq.0)then
	  H_mol=0
	  return
	endif
	lgn=log10(n)
	lgT=log10(T)
	x_n=(lgn+2)*64/13
	x_T=lgT*1024/7
	i_n=int(x_n)
	if(i_n.lt.0)then
	  H_mol=0
	  return
	endif
	i_T=int(x_T)
	if(i_T.lt.0)then
	  H_mol=0
	  return
	endif
	if(i_n.gt.63)i_n=63
	if(i_T.gt.1023)i_T=1023
	H0=Hy(i_T,i_n)
	i_n1=i_n
	if(i_n.lt.63)i_n1=i_n+1
	D_H_n=Hy(i_T,i_n1)-H0
	i_T1=i_T
	if(i_T.lt.1023)i_T1=i_T+1
	D_H_T=Hy(i_T1,i_n)-H0
	D_n=x_n-i_n
	D_T=x_T-i_T
	D_H=D_H_n*D_n+D_H_T*D_T
	H=H0+D_H
	H_mol=H
	return
      end
!
      subroutine Get_Table
	implicit real (a-z)
	common/H_H2_CO/ Lamb(0:1023,0:63),Hy(0:1023,0:63)
	open(1,file='H+H2+CO.bin'
     &	,status='old',form='unformatted')
	read(1)Lamb
	close(1)
	return
	end
!
      function heating(n_H,n_HI)
	implicit real (a-z)
	parameter(HR=3.8e-29)
	parameter(Hd=2.2e-28)
	parameter(zH=0.1)
	!print*,n_H,' cm^-3',n_HI,' cm^-3'
	Gamm_HR=HR*n_H
	Gamm_Hd=Hd*zH*n_H*n_HI
	Gamma= Gamm_HR + Gamm_Hd
	!print*,Gamma,' ergs cm^-3 sec^-1'
	!warning: positively returned value.
	!*** heating = Gamma
	heating = Gamma
	!*** heating = 4*Gamma
	!*** heating = 8*Gamma
	!*** heating = 32*Gamma
	!*** heating = 48*Gamma
	return
	end
