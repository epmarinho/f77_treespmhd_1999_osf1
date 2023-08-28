#ifndef _N_BODY
!///////////////////////////////////////////////////////////////////////
!         SPH AND TREE-NEAREST-NEIGHBORS-SEARCHING SUBROUTINES
!///////////////////////////////////////////////////////////////////////

#	ifndef _NEW_H_SEARCH
      subroutine Retry_convergence(i)
      INCLUDE 'SPH_common.h'
      itercounter = 0
      iterplus = iterplus + 1
      y = y * guessing
      ncount = n_fix/2
      if(iterplus.gt.ncount)then
          guessing = .5*(ran(iseed(i))+.998)
          if (y.lt.1.0e-4) then
            y = pi*1.0e5*guessing*exp(1.)*(1+pi/10)
          endif
          iterplus = 0
      endif
      end


      subroutine Set_permit
      INCLUDE 'SPH_common.h'
      do node=1,maxnode
        permit(node)=.true.
      end do
      end


      subroutine NNS_gather(i)
      INCLUDE 'SPH_common.h'
      neighb(i) = 0
      call Set_permit
      do node=1,maxnode
        if(permit(node))then
          call Search_gather_neighbors(i,2*h(i))
        else
          call Forbid(node)
        endif
      end do
      end


      subroutine Search_gather_neighbors(i,r_i)
      /*
        This subroutine is nothing more than a tree-walk procedure to perform
        a searching for particles that are included within a radius r_i from
        the particle i.
      */
      integer cell
      INCLUDE 'SPH_common.h'
      r_cell = arest(node) + r_i
      do cell=1,8
        if(np(node,cell).gt.0)then
!          if a unit cell, identify this one as a neighbor candidate:
          if (np(node,cell).eq.1) then
            if(spacing(cell,i,r_i,1).lt.0)neighb(i)=neighb(i)+1
          else
            if (spacing(cell,i,r_cell,0).gt.0) then
!            forbid tree descents for all next searchings:
              next = down(node,cell)
              permit(next) = .false.
            endif
          endif
        endif
      end do
      return
      end



      subroutine Neighboring
      INCLUDE 'SPH_common.h'
      !print*,'Neighboring: start'
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        neighb(i) = 0
        call Set_permit
        do node=1,maxnode
          if(permit(node))then
            call SRN(i,2*h(i))
          else
            call Forbid(node)
          endif
        end do
      end do
      !print*,'Neighboring: done'
      end



      subroutine SRN(i,r_i) /* Search for Relevant Neighbors */
      /*
      WARNING:

      This method for relevant neighbors searching has a 
      serious "bug". If we worry on preserving the neighbors
      within max{2hi,2hj} radius, it may occurs situations for
      which many particles, mainly from rarefied regions, are 
      included in the list. So litle, or nothing, can be 
      done to prevent this situation. For example, if there 
      are 128 neighbors with radii large enough to englobing 
      the particle i, even if these ones are out of the 2hi 
      radius, then all these neighbors will be added up to 
      the neighboring list.
      */
      integer cell
      INCLUDE 'SPH_common.h'
      r_cell = arest(node) + r_i
      do cell=1,8
        if (np(node,cell).gt.0) then
          if (np(node,cell).eq.1) then
            j = label(node,cell)
            r = amax(r_i,2*h(j))
            if (spacing(cell,i,r,1).lt.0) then
            if (neighb(i).lt.maxneighb) then
                neighb(i) = neighb(i) + 1
              neighb_list(i,neighb(i)) = j
              endif
            endif
          else
            if (spacing(cell,i,r_cell,0).gt.0) then
              next = down(node,cell)
              permit(next) = .false.
            endif
          endif
        endif
      end do
      end


      subroutine Search_lengths
      INCLUDE 'SPH_common.h'
c     print*,'Search_lengths: start'
      do i_=1,n_p(timebin)
        call Particle_length ( p_list(i_,timebin) )
      end do
c     print*,'Search_lengths: done'
      end


      subroutine Particle_length(i)
      INCLUDE 'SPH_common.h'
!      initialize control variables:
      y           = 100
      guessing    = .3427301914159
      iterplus    = 0
      itercounter = 0
      converging  = .false.
!      start iterations counter:
      iter_neighb = 0
      do while(.not.converging)
!        predict the smoothing-length:
        h(i)= h(i)*(1+y*((1+y)*n_fix/(n_fix+neighb(i)*y))**.3333333)/(1+y)
        call NNS_gather(i)
        call Ask_for_convergence(i)
        iter_neighb = iter_neighb + 1
      end do
      !print*,i,h(i),neighb(i)
      end


      subroutine Ask_for_convergence(i)
      INCLUDE 'SPH_common.h'
      itercounter = itercounter + 1
      ovc = (itercounter.gt.3)
      nonzero=(neighb(i).gt.1).and.(neighb(i).le.maxneighb)
      ok=(abs(neighb(i)-n_fix).le.(tol_nfix))
      converging = ok.and.nonzero
      if(ovc)then
        if(.not.converging) call Retry_convergence(i)
      endif
      end


      function spacing(cell,i,r,metric)
      integer cell
      INCLUDE 'SPH_common.h'
      d_i_c = 0
      if(metric.eq.1)then
            ! square-metric:
            do k=1,dim
              a = x_cell(node,cell,k) - xp(i,k)
              d_i_c = d_i_c + a*a
            end do
            spacing=d_i_c-r*r
      else
            ! metric of maximum:
            do k=1,dim
              a = abs( x_cell(node,cell,k) - xp(i,k) )
            d_i_c = amax(d_i_c,a)
            end do
            spacing=d_i_c-r
      endif
      end







#	else

!///////////////////////////////////////////////////////////////////////
!////////////// New version
!///////////////////////////////////////////////////////////////////////
!///////////////////////////////////////////////////////////////////////
!     sph and tree-nearest-neighbors-searching subroutines
!///////////////////////////////////////////////////////////////////////
      subroutine search_lengths
      include 'SPH_common.h'
      print*,'search_lengths: start'
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        call smoothing_length(i)
      end do
      print*,'search_lengths: done'
      end

      subroutine smoothing_length(i)
      logical repeat
      include 'SPH_common.h'
#ifdef _DEBUG
      print*,'calculating h for poarticle',i
#endif
      itercount=0
      repeat=.true.
      do while(repeat)
        do while(repeat)
          itercount=itercount+1
          call nns_gather(i,repeat) !-> neighb(i), neighb_list(i,...)
        end do
        if(neighb(i).gt.n_fix)then
          call sort_list(i)
        else
          if(neighb(i).lt.n_fix)then
            h(i)=
     &      0.50356*
     &      (1+(2*float(n_fix)/(n_fix+neighb(i)))**0.3333333)*h(i)
            repeat=.true.
          endif
          if(repeat.and.(itercount.gt.8))then
            h(i)=2.0273541*h(i)
            itercount=0
          endif
        endif
      end do
#ifdef _DEBUG
      print*,'yielding',h(i)
#endif
      end

      subroutine nns_gather(i,error)
      include 'SPH_common.h'
      logical error
      error=.false.
      neighb(i) = 0
      call set_permit
      do node=1,maxnode
        if(permit(node))then
          call search_gather_neighbors(i,2*h(i),error)
          if(error)then
            h(i)=0.5*h(i)
            return
          endif
        else
          call forbid
        endif
      end do
      end

      subroutine search_gather_neighbors(i,r_i,error)
      logical error
      integer cell
      include 'SPH_common.h'
      r_cell = arest(node) + r_i
      do cell=1,8
        if(np(node,cell).gt.0)then
!          if this is a unit cell, identify this one as a neighbor
!         candidate:
          if(np(node,cell).eq.1)then
            if(spacing(cell,i,r_i,1).lt.0)then
                neighb(i)=neighb(i)+1
              if(neighb(i).gt.maxneighb)then
!               there are too much neighbors for this input radius (r_i)
                neighb(i)=0
                error=.true.
                return
              endif
              j=label(node,cell)
              if(j.gt.n)then
                print*,'error in particle index'
                stop
              endif
              neighb_list(i,neighb(i)) = j
            endif
          else
            if(spacing(cell,i,r_cell,0).gt.0)then
!                forbid tree descents for internal searchings:
              next = down(node,cell)
              permit(next) = .false.
            endif
          endif
        endif
      end do
      error=.false.
      end

      subroutine sort_list(i)
        include 'SPH_common.h'
        do j_=1,neighb(i)
          j=neighb_list(i,j_)
          r2=0
          do k=1,3
            r2=r2+(xp(i,k)-xp(j,k))*(xp(i,k)-xp(j,k))
          end do
          !this is the sort-key:
          t(j_,1)=int(0.5*r2*neighb(i)/(h(i)*h(i)))
          !this is the information
          t(j_,2)=j
        end do
        nsafe=n
        n=neighb(i)
        call hsort
        n=nsafe
        do j_=1,n_fix
          neighb_list(i,j_)=t(j_,2)
        end do
        j=neighb_list(i,n_fix)
        neighb(i)=n_fix
        r2=0
        do k=1,3
           r2=r2+(xp(i,k)-xp(j,k))*(xp(i,k)-xp(j,k))
        end do
        h(i)=0.5*sqrt(r2)
      end

      function spacing(cell,i,r,metric)
      integer cell
      include 'SPH_common.h'
      d_i_c = 0
c     quadratical metric
      if(metric.eq.1)then
        do k=1,dim
          a = x_cell(node,cell,k) - xp(i,k)
          d_i_c = d_i_c + a*a
        end do
        spacing=d_i_c-r*r
      else
c     metric of the maximum:
        do k=1,dim
          a = abs( x_cell(node,cell,k) - xp(i,k) )
          d_i_c = amax(d_i_c,a)
        end do
        spacing=d_i_c-r
      endif
      end

      subroutine effective_neighbors(i,r_i) 
c      warning:
c      this method for relevant neighbors searching has a 
c      serious "bug". if we worry on preserving the neighbors
c      within max{2hi,2hj} radius, it may occurs situations for
c      which many particles, mainly from rarefied regions, are 
c      included in the list. so litle, or nothing, can be 
c      done to prevent this situation. for example, if there 
c      are 128 neighbors with radii large enough to englobing 
c      the particle i, even if these ones are out of the 2hi 
c      radius, then all these neighbors will be added up to 
c      the neighboring list.
      integer cell
      include 'SPH_common.h'
      r_cell = arest(node) + r_i
      do cell=1,8
        if(np(node,cell).gt.0)then
          if(np(node,cell).eq.1)then
c           an individual particle is now foccused:
            j = label(node,cell)
c           define the effective smoothing radius as the maximum
c           of the individual particle radius: 2*h(i) and 2*h(j)
            r = amax(r_i,2*h(j))
c           this choice avoids the underestimating of the kernel
c           contribution:
            if(spacing(cell,i,r,1).lt.0)then
              if(neighb(i).lt.maxneighb)then
                neighb(i) = neighb(i) + 1
                neighb_list(i,neighb(i)) = j
              endif
            endif
          else
            if(spacing(cell,i,r_cell,0).gt.0)then
              next = down(node,cell)
              permit(next) = .false.
            endif
          endif
        endif
      end do
      end


      subroutine neighboring
      include 'SPH_common.h'
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        neighb(i) = 0
        call set_permit
        do node=1,maxnode
          if(permit(node))then
            call effective_neighbors(i,2*h(i))
          else
            call forbid
          endif
        end do
      end do
      end

      subroutine set_permit
      include 'SPH_common.h'
      do node=1,maxnode
        permit(node)=.true.
      end do
      end

#	endif


#endif
