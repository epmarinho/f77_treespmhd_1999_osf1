#ifdef  _SPH_ADIABAT
#define _ENERGY_CHECK
#endif

      subroutine Settings
      INCLUDE 'SPH_common.h'
      iiseed=417253423
#ifndef _N_BODY
#  ifdef _2D
      n_fix = 19
#  else
      n_fix = 90
#  endif
#  if defined ( _SPH_ADIABAT )
      alpha = 0.5
      beta  = 1
      eta   = 0.1
      courant = .0625
#  else
      alpha = 0.5
      beta  = 1.0
      eta   = 0.1
      courant = 2.
#  endif
      tol_nfix = 1
#  if  ! (defined ( _SPH_ADIABAT ) || defined ( _TUBE ))
      tol_nfix = max(int( 4./sqrt(float(n_fix)) ),1)
#    ifndef _2D
      sc = n_fix**(.3333333)
#    else
      sc = sqrt(float(n_fix))
#    endif
#  endif
#else
      courant = .5
      W_p = 0
#endif
      end

      subroutine TimeSet
      INCLUDE 'SPH_common.h'
      dtm=dt(0)
      do j=1,ntbins
        dtm=.5*dtm
      end do
#ifdef _VERBOSE
      print*,'dtm =',dtm
#endif
      dtm_safe=dtm
      end

      subroutine DoSimulation
      INCLUDE 'SPH_common.h'
      external max, min, amax, amin
      common/general_counter/ icounter
      icounter=0
#ifdef _MOVIE_
      print*,'Movie simulation'
      if(global_time_level.eq.1)call Movie (global_time_level)
#endif
      print*,'Settings'
      call Settings
      print*,'Preparations'
      call Preparations
      print*,'TimeSet'
      call TimeSet

!     recover intervals-between-frame:
      is = global_time_level + 1
      do while(is.gt.ns)
        is = is - ns
      end do

!     get initial-frame-number:
      ifile = ini

!     recover remainning number of integrations:
      ninew = global_time_level
      do while (ninew.ge.ni)
        ninew = ninew - ni
      end do
      ninew = ni - ninew

#if !(defined(_TUBE)||defined(_SPH_ADIABAT)||defined(_N_BODY)||defined(_SPH_ISOTHERM)||defined(_GALAXY_2D))
      print*,'Get_Table'
      call Get_Table
#endif
#if defined(_ENERGY_CHECK) && ( defined(_SPH_ADIABAT) || defined(_TUBE) ) 
      dE = 0
#endif

!     perform the simulation now:
      do time_level=1,ninew

          global_time_level = global_time_level + 1

          rtime = global_time_level*dt(0)

#ifdef _VERBOSE
          write(*,'(a6,i5,a2,f15.7,a5)')
     &            'Step',time_level,' (',rtime,' [t])'
#endif

#if (defined(_SPH_STR)||defined(_SPH_STR_SN))
          str_growth = 0
#endif

          if(dt(0).eq.0.0)print*,'simula: warning: dt(0) is null'
#if defined (_ENERGY_CHECK) && ( defined(_SPH_ADIABAT) || defined(_TUBE) ) 
          E_old = E_tot ! E_tot is calculated in Integral_quantities
          call Integration_scheme (dE)
          if(time_level.gt.1) dE = amax (16*(E_old-E_tot)/(abs(E_old)),0)
#   if defined (_VERBOSE)
          print*,'energy correction: ',dE
#   endif
#else
          call Integration_scheme (0)
          if(dt(0).eq.0.0)then
            call finalize
            stop
          endif
#endif

#if (defined(_SPH_STR)||defined(_SPH_STR_SN))
          n_str = n_str + str_growth
#endif

!          integral parameters saving:
!          ---------------------------

!          extract prefix from output-file-names:

          l3        = index(percent,' ')

          if((global_time_level.eq.1).or.(Eo.eq.0))then
             open(2,file=relat,status='unknown',access='append')
             write(2,*)'file ',percent(1:l3)
             write(2,*)'Total energy = ',E_tot
             write(2,*)'Total gravitational energy = ',W_tot
             Eo=E_tot
             E_max=Eo
             E_min=E_max
             close(2)
          else
            E_max=amax(E_max,E_tot)
            E_min=amin(E_tot,E_min)
            ee = (E_max-E_min)/abs(Eo)
            open(2,file=percent,status='unknown')
            write(2,*)'energy has changed within ',ee*1.0e+02,' %'
            write(2,*)'last step =',global_time_level
            write(2,*)'physical time =',rtime
            write(2,*)'bottom timebin =',n_timebins
            write(2,*)'deepest timebin =',deepest
            write(2,*)'last time histogram:'
            write(2,'(3x,i3,i6,e19.10)')
     &            (ihist,n_p(ihist),dt(ihist),ihist=1,n_timebins)
            write(2,*)
            close(2)
          endif

          open(2,file=evolution,status='unknown',access='append')
          write(2,'(x,14e16.8)')
     &             rtime                 ! 1
     &            ,W_tot+U_tot+T_tot     ! 2
     &            ,U_tot                 ! 3
     &            ,T_tot                 ! 4
     &            ,W_tot                 ! 5
#ifndef _N_BODY
     &            ,W_tot+W_s+W_p+2*T_tot ! 6
     &            ,luminosity            ! 7
#else
     &            ,W_tot+W_s+2*T_tot     ! 6
     &            ,0.0e+00               ! 7
#endif
     &            ,(vmoment(j),j=1,dim)  ! 8..10
#ifndef _TUBE
     &            ,(angmom(j), j=1,dim)  ! 11..13
#endif
#ifdef _SPMHD
     &            ,E_mag                 ! 14
#endif
          close(2)

#if !( defined (_SPH_ADIABAT) || defined (_TUBE) || defined (_N_BODY) )
#   if (defined(_SPH_STR)||defined(_SPH_STR_SN))
          open(2,file=namestar,status='unknown',access='append')
          write(2,'(x,3e16.8)')
     &            rtime,                 ! 1
     &            n_str,                 ! 2
     &            str_growth/dt(0)       ! 3
          close(2)
#   endif
#endif

!          definitive output:
#ifdef _MOVIE_
          call Movie (global_time_level)
#endif
          if(is.eq.ns)then
            call finalize
            is=1
            ifile = ifile + 1
          else
            is=is+1
          endif

#ifndef _GALAXY_2D
!         security saving:
!         ------------------
          call writecore
#endif

      end do

      print*
      print*,'Simulation done! Have a nice paper!'
      print*

      end

