# 1 "SPH_io.F"



! 
!
!      I/O AND GENERAL SETTING UP SUBROUTINES:
!
! 
!
!      1) CORE GENERATOR/READER:
!
! 

      subroutine writecore
      INCLUDE 'SPH_common.h'




      open(3,file=namecore,status='unknown',form='unformatted')



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
     &             epsilon      !--> softening length
      write(3)
     &             E_max
      write(3) 
     &             Eo
      write(3) 
     &             E_min

      write(3)
     &             n_fix
      write(3)
     &             alpha
      write(3)
     &             beta
      write(3)
     &             eta

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








      write(3) 
     &            sc





      ! print*,'trying to close ',namecore
      close(3)
      ! print*,'close successful.'
      ! print*,'done!'
      ! print*
      end

      subroutine readcore
      INCLUDE 'SPH_common.h'
      character *32 lntxt
      print*
      print*,'readcore'
      print *
      print *,'   name of core-file (with no file type):'
      read(*,'(132a)')lntxt
      read(*,'(132a)')namecore
      print *
      l         = index(namecore(1:132),' ')
      namecore  = namecore(1:l-1) 
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
     &             epsilon      !--> softening length
      read(1)
     &             E_max
      read(1) 
     &             Eo
      read(1) 
     &             E_min

      read(1)
     &             n_fix
      read(1)
     &             alpha
      read(1)
     &             beta
      read(1)
     &             eta

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








      read(1) 
     &             sc





      close(1)
      epsilon2=epsilon*epsilon
      theta2=theta*theta
      print*,n,' particles was read'
      print*,'done!'
      print*
      end


! 
!
!      INPUT PARAMETERS:
!
! 


      subroutine readparameters
      INCLUDE 'SPH_common.h'
      esc      =char(27)
      cr         =char(13)
      cutl      =esc 
      home      =esc 
      cls       =esc 
      rev       =esc 
      flsh      =esc 
      hltd      =esc 
      norm      =esc 
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












      !print*,'enter the background pressure/density (computer units):'
      Pbg=.87025



      print *
      print *,'parameters ok!'
      theta2=theta*theta




!     find a convenient time-setp:
      if(dt(0).lt.0)then







        print*,'invalid softening length. aborting.'
        stop

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

     &   ,u(i)                               ! 8




      if(i.eq.n)go to 11
      go to 5
   10      n=i-1
   11      close(1)
      do i=1,n
        Tot_Mass = Tot_Mass + mass(i)
        !print*,i,mass(i),(xp(i,j),j=1,dim),(v(i,j),j=1,dim),u(i)
      end do
      print *
      print *,n,' phase-points were read'
      print *
      end























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
      if(ichoice.eq.1)then
        call readparameters
      else
        call readcore
      endif
!
!      extract prefix from output-file-names:
!
      l1        = index(nameout,']')+1
      l2        = index(nameout(l1:132),'.')-1
      l3        = index(name,' ')
      l4        = index(nameout,' ')
      l5        = index(namecore,' ')
      namecore  = nameout(1:l1+l2) 
      evolution = nameout(1:l1+l2) 
      relat     = nameout(1:l1+l2) 
      percent   = nameout(1:l1+l2) 



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
      print *,'output file sequence             : ',nameout(1:l4)
      print *,'number of next output            : ',ini

      print *,'number of neighbors              : ',n_fix
      print *,'alpha                            : ',alpha
      print *,'beta                             : ',beta
      print *,'eta                              : ',eta
      print *
      print *



          print *,'      dissipative simulation!'




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

      write(4,*)'number of neighbors             : ',n_fix
      write(4,*)'alpha                           : ',alpha
      write(4,*)'beta                            : ',beta
      write(4,*)'eta                             : ',eta







      write(4,*)'      non-adiabatic simulation'


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
! 


      


      subroutine Finalize
      character*132 ext



      character*132 pre2




      character chrd*4
      INCLUDE 'SPH_common.h'
      l1  = index(nameout,']')+1
      l2  = index(nameout(l1:132),'.')-1
      l   = l1+l2
      ext = nameout(l:132)
      call numstring(chrd,mchr,ifile)



      pre2 = nameout(1:l-1) 

      open(2,file=pre2,status='unknown')
      print *,'       ',pre2
























      do i=1,n
        if(gas(i))then















              write(2,'(14e16.8)')


     &         mass(i)            !  1
     &        ,(xp(i,j),j=1,dim)  !  2 - 4
     &        ,(v(i,j),j=1,dim)   !  5 - 7
     &        ,u(i)               !  8



     &        ,rho(i)             !  9 ! 12
     &        ,P(i)*rho(i)*rho(i) ! 10 ! 13
     &        ,-lambda(i)         ! 11 ! 14
     &        ,Temp(i)            ! 12 ! 15
     &        ,cs(i)              ! 13 ! 16













        endif
      end do



      close(2)

      end




      subroutine numstring(chrd,mchr,ifile)
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
        chrd = char(idig+48) 
        itmp = itmp/10
      end do
      end



      subroutine Movie(gtl)
      integer gtl
      character*132 ext,pre4
      character chrd*4
      !*** integer leg
      INCLUDE 'SPH_common.h'
      !*** leg=mod(gtl,4)
      !*** if(leg.ne.0)return
      !*** print*,'NOTE: The animation file has been saved with the addition'
      !*** print*,'      of a density column.'
      ! npoints=4096
      npoints=16384
      
      l1  = index(nameout,']')+1
      l2  = index(nameout(l1:132),'.')-1
      l   = l1+l2
      ext = nameout(l:132)
      call numstring(chrd,mchr,gtl)
      pre4 = nameout(1:l-1) 
      print*,'saving movie file:',pre4
      n_show=n/npoints
      if(n_show.eq.0)n_show=1
      open(4,file=pre4,status='unknown')
      do i=1,n,n_show

          if(gas(i))then
             write(4,'(5e12.4)')
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
            write(4,'(5e12.4)')
     &       (xp(i,j),j=1,dim)
*    &      ,(v(i,j),j=1,dim)
*    &      ,u(i)
*    &      ,P(i)*rho(i)*rho(i)
*    &      ,-lambda(i)
*    &      ,1e+4
*    &      ,1e+4
          endif






      end do
      close(4)
      end

