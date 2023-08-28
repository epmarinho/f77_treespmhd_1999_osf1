# 1 "SPH_grav.F"
# 305 "SPH_grav.F"

      subroutine extern_gravity_2d
      INCLUDE 'SPH_common.h'
      !print*,'extern_gravity_2d'

      ! physical units must be given as [l]=1 kpc, [t]=100 Myr
      ! unperturbed parameters:
      v1=33.77 != 330 km/s
      v2=86.98 != 850 km/s
      ra1=18.0 !kpc
      rb1=1.50 !kpc
      ra2=0.73 !kpc
      rb2=0.30 !kpc
      rcore=0.01!kpc
      rcore=rcore*rcore

      ! perturbation parameters:
      pitch=-(7.0/180.0)*pi
      pitch=1/tan(pitch)
      R0=7.5 !kpc
      Ampl=-2.612 ! [v^2]
      r=R0
      Omega_p=2.4291!(v1*exp(-r/ra1-rb1/r)+v2*exp(-r/ra2-rb2/r))/r
      print*,'Omega_p=',Omega_p
      nspirals=2

      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        !print*,'i=',i
        r2=0
        do j=1,dim
          r2=r2+xp(i,j)*xp(i,j)
        end do
        r=sqrt(r2+rcore)
        !print*,'r=',r

        !*** unperturbed potential ***
        vc=v1*exp(-r/ra1-rb1/r)+v2*exp(-r/ra2-rb2/r)
        !print*,'vc=',vc
        gc=-vc*vc/r
        !print*,'g=',gc


        !*** spiral perturbation ***
        angle=acos(xp(i,1)/r)
        if(y.lt.0)angle=-angle
        print*,'angle=',angle
        gs=Ampl



     &    *sin(
     &      nspirals*(
     &        pitch*log(r/R0)-angle
!#    ifndef 1
     &        +Omega_p*time(timebin)
!#    endif  1
     &      )
     &    )
     &    /r
        g(i,1)=gs*(pitch*xp(i,1)+xp(i,2))/r
        g(i,2)=gs*(pitch*xp(i,2)-xp(i,1))/r


        !*** finalization
        do j=1,dim

          g(i,j)=g(i,j)+gc*xp(i,j)/r




          !*** include centrifugal force (don't forget the Coriolis one!)
          g(i,j)=g(i,j)+Omega_p**2*xp(i,j)

        end do

      end do
      !print*,'extern_gravity_2d: done'
      write(*,'(a16,e15.7,a11,i2)')
     &'potential time=',time(timebin),' time-bin=',timebin
      end


      subroutine Coriolis
      INCLUDE 'SPH_common.h'
        !*** calculate the Coriolis force:
        Omega_p=2.4291
        do i_=1,n_p(timebin)
          i=p_list(i_,timebin)
          accel(i,1)=accel(i,1)+2*v(i,2)*Omega_p
          accel(i,2)=accel(i,2)-2*v(i,1)*Omega_p
        end do
        !*** end of Coriolis force:
      end



