      program spheres
        real Mt
	include "iso_g.h"
!
	character answer
	pi = 4*atan(1.)
	print *,'arquivo '
	read(*,'(a128)')nameout
	print *,'system total mass'
	read(*,*)Mt
	print *,'tidal radius'
	read(*,*)Rt
	print *,'core radius'
	read(*,*)Rcore
	print*,'adiabatic constant (average):'
	read(*,*)gamma
	print*,'Magnetic field components (initially homogeneous):'
	read(*,*)(B(j),j=1,dim)
	print *,'numero de particulas'
	read(*,*)n
	print *,'Gerando ',n,' particulas.'
	print *,'semente'
	read(*,*)iseed
	iseed=iseed+(mod(iseed,2)-1)
	print*,'iseed =',iseed
	jseed=1402233233*ran(iseed)
	jseed=jseed+(mod(jseed,2)-1)
	print*,'jseed =',jseed
	kseed=2093233413*ran(jseed)
	kseed=kseed+(mod(kseed,2)-1)
	print*,'kseed =',kseed
	call spher_2(Mt,Rt,Rcore)
	call coordinates
	open(1,file=nameout,status='unknown',access='append')
	mass = Mt/float(n)
	do i=1,n
	  write(1,'(11e16.8)')
     &	  mass,
     &	  (x(i,j),j=1,dim),(v(i,j),j=1,dim),u(i),(B(j),j=1,dim)
	end do
	close(1)
	l2        = index(nameout,' ') - 1
	open(1,file=nameout(1:l2)//'_head',status='unknown')
	write(1,*)'arquivo ',nameout(1:l2)
	close(1)
      end


      subroutine spher_2(Mt,Rt,Rcore)
c	1/r sphere
        real Mt
	include 'iso_g.h'
	Rt_2=Rt*Rt
	do i=1,n
	  mass = ran(iseed)
	  r(i) = Rt*sqrt(mass)
	  !u(i) = Mt/(2*pi*Rt_2*r(i))
	  u(i) = 1/(gamma-1)*(Mt/Rt_2)*r(i)*abs(log(Rt/(r(i)+Rcore)))
	end do
	rho_m=Mt/(2*pi*Rt_2*Rt)
      end

      subroutine coordinates
	include "iso_g.h"
	pi = 4*atan(1.)
	do i=1,n
	  chi = ran(jseed)
	  !write(*,'(a4,i12,a2,e13.4)')'ran(',jseed,')=',chi
	  ! azimuthal angle is straightly calculated:
	  phi = 2*pi*ran(kseed)
	  ! cos(theta):
	  costheta = 1-2*chi
	  sintheta = 2*sqrt(chi*(1-chi))
	  ! xy-projection:
	  x(i,1)   = r(i)*sintheta
	  ! z-projection:
	  x(i,3)   = r(i)*costheta
	  ! y-projection:
	  x(i,2)   = x(i,1)*cos(phi)
	  ! x-projection:
	  x(i,1)   = x(i,1)*sin(phi)
	end do
      end
