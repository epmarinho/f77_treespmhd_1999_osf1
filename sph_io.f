!///////////////////////////////////////////////////////////////////////
!
!      I/O AND GENERAL SETTING UP SUBROUTINES:
!
#ifndef _GALAXY_2D
!///////////////////////////////////////////////////////////////////////
!
!      1) CORE GENERATOR/READER:
!
!///////////////////////////////////////////////////////////////////////

      subroutine writecore
      INCLUDE 'SPH_common.h'
#ifdef _VERBOSE
      print*
      print*,'writing over ',namecore
#endif
      open(3,file=namecore,status='unknown',form='unformatted')
#ifdef _VERBOSE
      print*,'successfully opening ',namecore
#endif
      write(3)
     &             name
      write(3)
     &             nameout
      write(3)
     &             n            !--> number of computer-particles 
      write(3)
     &             ifile            !--> next frame number
      write(3)
     &             ni            !--> number of root-integrations
      write(3)
     &             ns            !--> number of steps by frame
      write(3)
     &             Tot_Mass
      write(3)
     &             theta            !--> apperture parameter
      write(3)
     &             epsilon
      write(3)
     &             E_max
      write(3) 
     &             Eo
      write(3) 
     &             E_min
#ifndef _N_BODY
      write(3)
     &             n_fix
      write(3)
     &             alpha
      write(3)
     &             beta
      write(3)
     &             eta
#endif
      write(3)
     &             global_time_level
      write(3)
     &             ntbins2
      write(3)
     &             n_timebins
      write(3)
     &             deepest
      write(3)
     &             dtm
      write(3)
     &             dt(0)
      write(3)
     &            dt_old
      write(3)
     &            mass
      write(3)
     &            xp
      write(3)
     &            v
#ifndef _N_BODY
      write(3)
     &            accel
      write(3)
     &            u
      write(3)
     &            u_dot
      write(3)
     &            vmonag
      write(3)
     &            gas
#if !(defined(_SPH_ADIABAT))&&(defined(_SPH_STR)||defined(_SPH_STR_SN))
      write(3)
     &            star_mass
      write(3)
     &            t_SN
      write(3) 
     &            n_str
#endif
      write(3) 
     &            sc
#ifdef _SPMHD
      write(3)
     &            B
#endif
#endif /* _N_BODY */
      ! print*,'trying to close ',namecore
      close(3)
      ! print*,'close successful.'
      ! print*,'done!'
      ! print*
      end

      subroutine readcore
      INCLUDE 'SPH_common.h'
      character *32 lntxt
      !print*
      !print*,'readcore'
      !print *
      !print *,'   name of core-file (with no file type):'
      read(*,'(132a)')lntxt
      read(*,'(132a)')namecore
      !print *
      l         = index(namecore(1:132),' ')
      namecore  = namecore(1:l-1)//'.core'
      open(1,file=namecore,status='old',form='unformatted')
      read(1)
     &             name
      read(1)
     &             nameout
      read(1)
     &             n            !--> number of particles
      read(1)
     &             ini            !--> initial frame-number
      read(1)
     &             ni            !--> number of root-integrations
      read(1)
     &             ns            !--> number of steps by frame
      read(1)
     &             Tot_Mass
      read(1)
     &             theta            !--> apperture parameter
      read(1)
     &             epsilon
      read(1)
     &             E_max
      read(1) 
     &             Eo
      read(1) 
     &             E_min
#ifndef _N_BODY
      read(1)
     &             n_fix
      read(1)
     &             alpha
      read(1)
     &             beta
      read(1)
     &             eta
#endif
      read(1)
     &             global_time_level
      read(1)
     &             ntbins2
      read(1)
     &             n_timebins
      read(1)
     &             deepest
      read(1)
     &             dtm
      read(1)
     &             dt(0)
      read(1)
     &            dt_old
      read(1)
     &            mass
      read(1)
     &            xp
      read(1)
     &            v
#ifndef _N_BODY
      read(1)
     &            accel
      read(1)
     &            u
      read(1)
     &            u_dot
      read(1)
     &            vmonag
      read(1)
     &            gas
#if !(defined(_SPH_ADIABAT))&&(defined(_SPH_STR)||defined(_SPH_STR_SN))
      read(1)
     &            star_mass
      read(1)
     &            t_SN
      read(1) 
     &             n_str
#endif
      read(1) 
     &             sc
#ifdef _SPMHD
      read(1)
     &             B
#endif
#endif
      close(1)
      theta2=theta*theta
      eps2_max=epsilon*epsilon
      eps2_min=eps2_max*float(n)**(-.3333333)
      print*,n,' particles was read'
      print*,'done!'
      print*
      end
#endif _GALAXY_2D


!///////////////////////////////////////////////////////////////////////
!
!      INPUT PARAMETERS:
!
!///////////////////////////////////////////////////////////////////////


      subroutine readparameters
      INCLUDE 'SPH_common.h'
      esc      =char(27)
      cr         =char(13)
      cutl      =esc//'[k'
      home      =esc//'[2h'
      cls       =esc//'[2j'
      rev       =esc//'[7m'
      flsh      =esc//'[5m'
      hltd      =esc//'[1m'
      norm      =esc//'[0m'
      deepest = 0
      print *
      print *,'system set up'
      print *
      print *,' file for initial conditions:'
      read(*,'(a132)')lntxt
      read(*,'(a132)')name
      print *,' enter N (0 for the maximum of the file, noting that'
     &       ,' the default',nn,' is the maximum):'
      read(*,'(a132)')lntxt
      read(*,*)n
      if((n.eq.0).or.(n.gt.nn))n=nn
      print *
      print *,'   generic output frame-name (with file type):'
      read(*,'(a132)')lntxt
      read(*,'(a132)')nameout
      print *
      ini=-1
      do while(ini.lt.0)
        print *,'   initial frame number:'
        read(*,'(a132)')lntxt
        read(*,*)ini
        if(ini.lt.1)then
          print*,'initial file number must be an unsigned integer'
        endif
      end do
      print *,'   up to final number of integrations:'
      read(*,'(a132)')lntxt
      read(*,*)ni
      print *,'   number of steps between frames:'
      read(*,'(a132)')lntxt
      read(*,*)ns
      print *
      print *,' last step number (zero to begin a simulation):'
      read(*,'(a132)')lntxt
      read(*,*)global_time_level 
      print *
      call Read_phase
#if ! defined (_NON_SELFGRAVITATING)
      print *
      print *,'enter tree-gravity parameters:'
      print *
      print *,'   aperture parameter:'
      read(*,'(a132)')lntxt
      read(*,*)theta
      if(theta.gt.1)then
        print *,'?: error: invalid theta; must be less than 1.'
        stop
      endif
#  ifdef _CONSTANT_SOFTENING
      print *,'   softening length (negative for default):'
      read(*,'(a132)')lntxt
      read(*,*)epsilon
#  endif
#endif
      print *
      print *,'integration parameters:'
      print *
      print *,'   root time-step (negative for default):'
      read(*,'(a132)')lntxt
      read(*,*)dt(0)
      print *,'   maximum number of time-bins (availb.=',ntbins,'):'
      read(*,'(a132)')lntxt
      read(*,*)ntbins2
      if(ntbins2.gt.ntbins)then
        ntbins2=ntbins
      else
        if(ntbins2.lt.1)then
          print *,'?: error: invalid number of time-bins !'
          stop
        endif
      endif
      print *
#ifndef _N_BODY
# ifdef _SPMHD
#  ifdef _BACKGROUND_FIELD_
      print*,'enter the background (homogeneous) field components ',
     &'(computer units):'
      read(*,'(a132)')lntxt
      read(*,*)(Bbg(i),i=1,dim)
#  endif
# endif
# ifdef _BACKGROUND_PRESSURE_
#  ifndef _GALAXY_2D
      print*,'enter the background pressure (computer units):'
      read(*,'(a132)')lntxt
      read(*,*)Pbg
#  else
      !print*,'enter the background pressure/density (computer units):'
      Pbg=.87025
#  endif
# endif
#endif
      print *
      print *,'parameters ok!'
      theta2=theta*theta
# ifndef _NON_SELFGRAVITATING
#   ifdef _CONSTANT_SOFTENING
      if(epsilon.lt.0) call calc_epsilon
#   else
      call calc_epsilon
#   endif
# endif
!     find a convenient time-setp:
      if(dt(0).lt.0)then
# ifndef _NON_SELFGRAVITATING
        dtm = epsilon*sqrt(epsilon/Tot_mass)*float(N)**(1./6)
        dt(0)    = 1024.
        do while(dt(0).gt.dtm)
          dt(0) = .5 * dt(0)
        end do
# else
        print*,'invalid time-step. aborting.'
        stop
# endif
      endif
!     determine the minimum admissible time-step:
*     dtm=2*dt(0)
*     do j=1,ntbins2
*       dtm=.5*dtm
*     end do
!     alternatively, you can make:
      dtm=2*dt(0)
      do j=1,ntbins
        dtm=.5*dtm
      end do
!     or yet (with your own risk!):
      dtm=2*dt(0)
      do j=1,11
        dtm=.5*dtm
      end do
      print*,'dtm =',dtm
      end

      subroutine Read_phase
      INCLUDE 'SPH_common.h'
      print*,'trying to open file:',name
      open(1,file=name,status='old')
      i=0
    5 i=i+1
      read(1,*,end=10)
     &    mass(i)                            ! 1
     &   ,(xp(i,j),j=1,dim),(v(i,j),j=1,dim) ! 2 - 7
#ifndef _N_BODY
     &   ,u(i)                               ! 8
#  ifdef _SPMHD
     &   ,(B(i,j),j=1,dim)                   ! 9 - 11
#  endif
#endif
      if(i.eq.n)go to 11
      go to 5
   10      n=i-1
   11      close(1)
      do i=1,n
        Tot_Mass = Tot_Mass + mass(i)
        !print*,i,mass(i),(xp(i,j),j=1,dim),(v(i,j),j=1,dim),u(i)
      end do
#     if defined(_ROTATING_FRAME)&&defined(_GALAXY_2D)
      Omega_p=2.4291
      do i=1,n
        v(i,1)=v(i,1)+Omega_p*xp(i,2)
        v(i,2)=v(i,2)-Omega_p*xp(i,1)
      end do
#     endif
      print *
      print *,n,' phase-points were read'
      print *
      end

# ifndef _NON_SELFGRAVITATING
      subroutine calc_epsilon
      logical cont
      INCLUDE 'SPH_common.h'
      print*,'calc_epsilon: start'
      gamma_ = -.25*Tot_Mass*Tot_Mass*float(n)**(-.3333333)
      cont = .true.
      epsilon = 0
      call maketree
      call treepotential
      epsilon  = gamma_/W_tot
      eps2_max = epsilon*epsilon
      eps2_min = eps2_max*float(n)**(-.3333333)
      do while(cont)
        epsilon_old = epsilon
        w_old = W_tot
        call treepotential
        epsilon  = gamma_/(W_tot+w_old)
        cont = abs(epsilon-epsilon_old).gt.(epsilon/32.)
        eps2_max = epsilon*epsilon
        eps2_min = eps2_max*float(n)**(-.3333333)
        print*,'W_tot = ',W_tot,' W_old = ',W_old
      end do
      print*,'epsilon_max^2 = ',eps2_max,' epsilon_min^2 = ',eps2_min
      print*,'calc_epsilon: done'
      end
# endif
!
!
!
!
!
!
!
!
      subroutine Preparations
      logical ciclica
      INCLUDE 'SPH_common.h'
      print *,'   input (1) for initializing or any other number to'
     &      ,' continue another simulation:'
      read(*,'(a132)')lntxt
      read(*,*)ichoice
      print *
!
      first=.true.
      not_yet=.true.
      starting=.true.
#ifndef _GALAXY_2D
      if(ichoice.eq.1)then
        call readparameters
      else
        call readcore
      endif
#else
      call readparameters
#endif
!
!      extract prefix from output-file-names:
!
      l1        = index(nameout,']')+1
      l2        = index(nameout(l1:132),'.')-1
      l3        = index(name,' ')
      l4        = index(nameout,' ')
      l5        = index(namecore,' ')
      namecore  = nameout(1:l1+l2)//'core'
      evolution = nameout(1:l1+l2)//'evol'
      relat     = nameout(1:l1+l2)//'rel'
      percent   = nameout(1:l1+l2)//'prcnt'
#if !(defined(_N_BODY)||defined(_SPH_ADIABAT))&&(defined(_SPH_STR)||defined(_SPH_STR_SN))
      namestar  = nameout(1:l1+l2)//'stars'
#endif
!
!      parameters confirmation
!
      print *
      print *
      print *,'number of particles              : ',n
      print *,'system mass                      : ',Tot_Mass
      print *,'aperture parameter               : ',theta
      print *,'softening length                 : ',epsilon
      print *,'time-step                        : ',dt(0)
      print *,'number of time-bins              : ',ntbins2
      print *,'number of root-integrations      : ',ni
      print *,'number of the last integration   : ',global_time_level
      print *,'number of steps between output   : ',ns
      print *,'initial conditions file-name     : ',name(1:l3)
      print *,'base name output sequence        : ',nameout(1:l4)
      print *,'number of next output            : ',ini
#ifndef _N_BODY
      print *,'number of neighbors              : ',n_fix
      print *,'alpha                            : ',alpha
      print *,'beta                             : ',beta
      print *,'eta                              : ',eta
      print *
      print *
#ifdef _SPH_ADIABAT
          print *,'      adiabatic simulation!'
#else
          print *,'      dissipative simulation!'
#endif
#else
      print *,'      pure gravitational simulation!'
#endif
      print *

      open(4,file=relat,status='unknown')
      write(4,*)'number of particles              : ',n
      write(4,*)'system mass                      : ',Tot_Mass
      write(4,*)'aperture parameter               : ',theta
      write(4,*)'softening length                 : ',epsilon
      write(4,*)'time-step                        : ',dt(0)
      write(4,*)'number of time-bins              : ',ntbins2
      write(4,*)'number of root-integrations      : ',ni
      write(4,*)'number of the last integration   : ',global_time_level
      write(4,*)'number of steps between output   : ',ns
      write(4,*)'initial conditions file-name     : ',name(1:l3)
      write(4,*)'output file sequence             : ',nameout(1:l4)
      write(4,*)'number of next output            : ',ini
#ifndef _N_BODY
      write(4,*)'number of neighbors             : ',n_fix
      write(4,*)'alpha                           : ',alpha
      write(4,*)'beta                            : ',beta
      write(4,*)'eta                             : ',eta
#  if defined(_SPMHD)
#    if defined(_BACKGROUND_FIELD_)
      write(4,*)'background magnetic field       : ',(Bbg(i),i=1,dim)
#    endif
#  endif
#  if defined(_BACKGROUND_PRESSURE_)
      write(4,*)'background pressure             : ',Pbg
#  endif
#ifdef _SPH_ADIABAT
      write(4,*)'      adiabatic simulation'
#else
      write(4,*)'      non-adiabatic simulation'
#endif
#endif
      close(4)
!
!
!   determination of levi-civita's unity third-order tensor.
!
      do i=1,dim
        do j=1,dim
          do k=1,dim
            if((i*j*k).eq.6)then
            l=i+3*(j+3*k)
                  ciclica=(l.eq.20).or.(l.eq.24).or.(l.eq.34)
            if(ciclica)then
              levi_civita(i,j,k) =  1
            else
              levi_civita(i,j,k) = -1
            endif
            else
                  levi_civita(i,j,k)   =  0
            endif
          end do
        end do
      end do
!
      end
!    
!///////////////////////////////////////////////////////////////////////


      /*            TIME-FINALIZATION:             */


      subroutine Finalize
      character*132 ext
#ifdef _N_BODY
      character*132 pre1
#else
      character*132 pre2
#if defined(_SPH_STR)||defined(_SPH_STR_SN)
      character*132 pre3
#endif
#endif
      character chrd*4
      INCLUDE 'SPH_common.h'
      l1  = index(nameout,']')+1
      l2  = index(nameout(l1:132),'.')-1
      l   = l1+l2
      ext = nameout(l:132)
      call numstring(chrd,mchr,ifile)

#if !defined(_N_BODY)

      pre2 = nameout(1:l-1)//'_gas'//chrd(1:mchr)//ext

      open(2,file=pre2,status='unknown')
      print *,'       ',pre2

#	if defined(_SPH_STR)||defined(_SPH_STR_SN)
      pre3 = nameout(1:l-1)//'_strs'//chrd(1:mchr)//ext
      open(3,file=pre3,status='unknown')
      print *,'       ',pre3
#	endif

#else
      pre1 = nameout(1:l-1)//chrd(1:mchr)//ext
      open(1,file=pre1,status='unknown')
      print *
      print *,'saving ',pre1
      do i=1,n
          write(1,'(7e16.8)')
     &     mass(i)            ! 1
     &    ,(xp(i,j),j=1,dim)      ! 2 ... 4
     &    ,(v(i,j),j=1,dim)      ! 5 ... 7
      end do
      close(1)

#endif

#if !defined(_N_BODY)

      do i=1,n
        if(gas(i))then

#      if !( defined (_H_2) || defined (_H_1) )

          rho_CGS=rho(i)/u_dens
          call y_weight(rho_CGS*XX,Temp(i),HI)! --> HI
          /* w_m = w_mol(XX,YY,HI) */

#      endif

#      if defined (_SPMHD)

              write(2,'(17e16.8)')

#      else

              write(2,'(14e16.8)')

#      endif
     &         mass(i)            !  1
     &        ,(xp(i,j),j=1,dim)  !  2 - 4
     &        ,(v(i,j),j=1,dim)   !  5 - 7
     &        ,u(i)               !  8
#      if defined (_SPMHD)
     &        ,(B(i,j),j=1,dim)   !  9 - 11
#      endif
     &        ,rho(i)             !  9 ! 12
     &        ,P(i)*rho(i)*rho(i) ! 10 ! 13
     &        ,-lambda(i)         ! 11 ! 14
     &        ,Temp(i)            ! 12 ! 15
     &        ,cs(i)              ! 13 ! 16
#      if !( defined (_H_2) || defined (_H_1) )
     &        ,HI                 ! 14 ! 17
#      endif
#      if defined(_SPH_STR) || defined(_SPH_STR_SN)
        else
              write(3,'(9e16.8,1x,i6)')
     &         mass(i)            ! 1
     &        ,(xp(i,j),j=1,dim)  ! 2 - 4
     &        ,(v(i,j),j=1,dim)   ! 5 - 7
     &        ,star_mass(i)/solar ! 8
     &        ,-lambda(i)/u_cool  ! 9
     &        ,i                  ! 10
#      endif
        endif
      end do
#      if defined(_SPH_STR) || defined(_SPH_STR_SN)
      close(3)
#      endif
      close(2)
#endif
      end




      subroutine numstring (chrd,mchr,ifile)
      character chrd*4
      if(ifile.gt.0)then
            mchr=log10(float(ifile))+1
      else
            mchr=1
      endif
      itmp=ifile
      chrd=char(0)
      do ichr=1,mchr
        idig = mod(itmp,10)
        chrd = char(idig+48)//chrd
        itmp = itmp/10
      end do
      end


#ifdef _MOVIE_
      subroutine Movie(gtl)
      integer gtl
      character*132 ext,pre4
      character chrd*4
      integer leg,iskip
      parameter(iskip=1)
      INCLUDE 'SPH_common.h'
      leg=mod(gtl,iskip)
      if(leg.ne.0)return
      leg=gtl/iskip
      ! npoints=4096
      npoints=16384
      /* npoints=n */
      l1  = index(nameout,']')+1
      l2  = index(nameout(l1:132),'.')-1
      l   = l1+l2
      ext = nameout(l:132)
      call numstring(chrd,mchr,leg)
      pre4 = nameout(1:l-1)//'_movie'//chrd(1:mchr)//ext
      print*,'saving movie file:',pre4
      n_show=n/npoints
      if(n_show.eq.0)n_show=1
      open(4,file=pre4,status='unknown')
      do i=1,n,n_show
#  ifndef _N_BODY
          if(gas(i))then
             write(4,'(3e12.4)')
     &        (xp(i,j),j=1,dim)
*    &       ,(v(i,j),j=1,dim)
*    &       ,u(i)
*    &       ,P(i)*rho(i)*rho(i)
*    &       ,-lambda(i)
*    &       ,Temp(i)
*    &       ,rho(i)
          endif
      end do
      do i=1,n,n_show
          if(.not.gas(i))then
            write(4,'(3e12.4)')
     &       (xp(i,j),j=1,dim)
*    &      ,(v(i,j),j=1,dim)
*    &      ,u(i)
*    &      ,P(i)*rho(i)*rho(i)
*    &      ,-lambda(i)
*    &      ,1e+4
*    &      ,1e+4
          endif
#      else
          write(4,'(3e11.3)')
     &           (xp(i,j),j=1,dim)
#      endif
      end do
      close(4)
      end
#endif
