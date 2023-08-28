      program perform_cooling_table
!
	implicit real*8 (a-z)
	integer l,n_rho,n_T,mchr,n_grid,flag
	character*32 nameout,ext,chrd,file_name
	parameter(T0=5.2490d4,n0=1.27022d+24,n_grid=1023)
!
	nameout='H+H2+CO.dat'
	!nameout='cooling.dat'
	open (1,file=nameout,status='unknown')
	!do n_rho=-2,8
	do n_rho=0,63

	  !n=(1d1)**n_rho
	  n=(1d1)**(13d0*n_rho/63-2)

	  !------------------------------------------
	  !l  = index(nameout(1:32),'.')
	  !ext = nameout(l:32)
	  !call numstring(chrd,mchr,n_rho+2)
	  !file_name=nameout(1:l-1)//chrd(1:mchr)//ext
	  !print*,'opening file ',file_name
	  !open (1,file=file_name,status='unknown')
	  !------------------------------------------

	  flag=1

	  do n_T=0,n_grid

	    T=10**( n_T*7d0/n_grid )
	    H = y(n,T,H2)

	    n_H2=(H2/(1+H))*n
	    n_H=n-n_H2

	  !-----------------------------------------
**	    if(n_H.gt.0)
**     &write(1,'(4e19.10)')dlog10(T),dlog10(n_H)
*	    if(n_H2.gt.0)
*     &	write(1,'(4e19.10)')dlog10(T),dlog10(n_H2)
*	write(1,'(4e19.10)')dlog10(n),dlog10(T),H,H2
	  !-----------------------------------------

	    Lamb=0

	    call cooling_H_H2(T,n_H,n_H2,L_H,L_H2)
	    Lamb=n_H2*(n_H2*L_H2+n_H*L_H)

	    call cooling_H_He(T,n_H,L_He)
	    Lamb = Lamb + L_He

	    Lamb = Lamb + (8.5d-5*n)*n*CO(T,n)

!	    if(n.gt.0)then

	      !Lamb = Lamb/n/n

!	      if(.not.(Lamb.gt.0))Lamb=1.0d-35
!      write(1,'(4e17.7)')dlog10(n),dlog10(T),dlog10(Lamb),H2
	      if(.not.(Lamb.gt.0))Lamb=1.0d-35
      write(1,'(4e17.7)')dlog10(n),dlog10(T),dlog10(Lamb)
!	      if(Lamb.gt.0)
!     &write(1,'(4e17.7)')dlog10(n),dlog10(T),dlog10(Lamb),H2
!			     1         2          3
!	    endif
	  end do
	  !close(1)
	end do
	close(1)
      end
!
      function y(n,T,H2)
	implicit real*8 (a-z)
	parameter(T0=5.249d+04,n0=6.3511d+23)
	T1=T0/T
	A=n0*dexp(-T1)
	A=A/n
	if(A.gt.1)then
	  y_aux=1/dsqrt(1+1/A)
	else
	  y_aux=dsqrt(A/(A+1))
	endif
	if((y_aux.gt.1).or.(y_aux.lt.0))then
	      print*,'y: error: y must be 0<= y <= 1'
	      print*,'y =',y_aux
	      print*,'A =',A
	      print*,'T =',T
	      stop
	endif
	H2=1-y_aux
	if(y_aux.eq.1)then
	  H2=.5/A
	endif
	y=y_aux
	return
      end
!
	!INCLUDE 'cooling.f'
	!INCLUDE 'CO.f'
	!INCLUDE 'numstring.f'
	!INCLUDE 'H-H2-HM.f'
