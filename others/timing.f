# 1 "SPH_timing.F"


      subroutine Start_list
      INCLUDE 'SPH_common.h'
      timebin=1
      dt(timebin)=dt(timebin-1)
      call Reset_time_list
      end

      subroutine multi_leapfrog (dE) !(FORTRAN IV-like approach)
      ! Here recursive calls are emulated by iterative
      ! procedures.
      INCLUDE 'SPH_common.h'
      timebin=1
 1110 continue!--> entry point
      dt_elapsed(timebin) = 0
 1111 continue!--> while (dt_elapsed(timebin).ne.dt(timebin-1)) do
      if(dt_elapsed(timebin).eq.dt(timebin-1)) goto 1114
          
          call Predictor_Corrector (dE)
          if(timebin.eq.n_timebins) goto 1113
              timebin = timebin + 1
              goto 1110!--> call Multiple_leapfrog:
 1112         continue !--> return point
 1113     continue!--> end if (.not.(timebin.eq.n_timebins))
          
          dt_elapsed(timebin) = dt_elapsed(timebin) + dt(timebin)
          time(timebin) = time(timebin) + dt(timebin)
          !print*,'dt_elapsed=',dt_elapsed(timebin),' timebin=',timebin
      goto 1111
 1114 continue!--> end while
      timebin = timebin - 1
      if(timebin.eq.0)return
      goto 1112!--> recursive return:
      end

      subroutine recover_dtm
      INCLUDE 'SPH_common.h'
      do while(dtm.lt.dtm_safe)
        courant=courant*2
        dtm=dtm*2
      end do
      end

      subroutine Reset_time_list
      integer i_p,top
      INCLUDE 'SPH_common.h'
      top    = timebin
      do while(n_timebins.gt.top)
        do i_p=1,n_p(n_timebins)
          p_list(i_p+n_p(top),top)=p_list(i_p,n_timebins)
        end do
        n_p(top)=n_p(top)+n_p(n_timebins)
        n_p(n_timebins) = 0
        n_timebins=n_timebins-1
      end do
      call Split_in_time_bins
      if(n_timebins.gt.deepest)deepest=n_timebins
      end

      subroutine Split_in_time_bins
      integer tbin,i_p,n_l,iparticle,top,bottom
      INCLUDE 'SPH_common.h'
      n_l = n_p(n_timebins) 
      n_p(n_timebins) = 0   
      bottom = n_timebins
      top    = bottom
      do i_p = 1, n_l
          iparticle = p_list(i_p, top)
          tbin      = top 
          dta       = dt(tbin)
          dtc       = dt_Courant(iparticle)
          do while ( (dta.gt.dtc) .and. (dta.gt.dtm) )
            dta  = .5 * dta
            tbin = tbin + 1
          end do
        tbin     = min(tbin,ntbins2)
          dt(tbin) = dta
          call put_in_daughter(iparticle,tbin,bottom)
      end do
      call Strip_p_list (bottom)
      end

      subroutine put_in_daughter(i,tbin,bottom)
      integer tbin,bottom,i
      INCLUDE 'SPH_common.h'
      n_p(tbin)=n_p(tbin)+1
      p_list(n_p(tbin),tbin) = i
      if(tbin.gt.bottom)bottom=tbin
      end

      subroutine Strip_p_list(bottom)
!            bottom is the bottom time bin.
!            up to now, n_timebins has not changed!
      integer tbin,bottom,pail,j
      INCLUDE 'SPH_common.h'
!      Strip list family by lines:
      if(bottom.gt.n_timebins)then
        pail=n_timebins-1
        do tbin=n_timebins,bottom
          if(n_p(tbin).gt.0)then
            pail=pail+1
!            copy tbin-row to pail:
            do j=1,n_p(tbin)
              p_list(j,pail)=p_list(j,tbin)
            end do
            dt(pail)=dt(tbin)
            n_p(pail)=n_p(tbin)
          endif
        end do
!        now, n_timebins is changed:
        n_timebins=pail
!        clear non-used quantities:
        if(pail.lt.ntbins2)then
          do tbin=pail+1,ntbins2
            dt(tbin)=0
            n_p(tbin)=0
          end do
        endif
      endif
      end

      function dt_Courant(i)
      INCLUDE 'SPH_common.h'
      vmax = 0
      gmax = 0
      smin=s(i)
      smin=courant*smin
      do j=1,dim
        vx=abs(v(i,j))
        gx=abs(accel(i,j))
        vmax=amax(vx,vmax)
        gmax=amax(gx,gmax)
      end do
      ttc = 64*dt(0)
      control=1
      do while ((control.gt.0).and.(ttc.gt.dtm))
        ttc=.5*ttc
        tau=.5*(ttc+dt_old(i))
        dtau=.25*(ttc-dt_old(i))
        control = (vmax+gmax*dtau)*tau - smin
      end do
      ttc=amax(ttc,1.0e-8)
      if(gas(i))then
        ttc=1/ttc
        tth=0
        if(dt_old(i).gt.0) tth=8*abs(h_dot(i))/dt_old(i)
        dtCourant = 
     &            (
     &              ttc +
     &              tth +
     &              (cs(i)+sqrt(vmonag(i)))/h(i)
     &            )
        dt_Courant = courant / dtCourant
      else
        dt_Courant = ttc
      endif
      end

      subroutine save_old_timesteps
      integer il,i
      INCLUDE 'SPH_common.h'
      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        dt_old(i) = dt(timebin)
      end do
      end
