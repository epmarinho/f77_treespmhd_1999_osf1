!///////////////////////////////////////////////////////////////////////
!     sph and tree-nearest-neighbors-searching subroutines
!///////////////////////////////////////////////////////////////////////
      subroutine search_lengths
            include 'sph_common.h'
          print*,'search_lengths'
            do i_=1,n_p(timebin)
                  i=p_list(i_,timebin)
                  call smoothing_length(i)
            end do
      end

      subroutine smoothing_length(i)
      logical repeat
      include 'sph_common.h'
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
      end

      subroutine search_gather_neighbors(i,r_i,error)
      logical error
      integer cell
      include 'sph_common.h'
      r_cell = arest(node) + r_i
      do cell=1,8
        if(np(node,cell).gt.0)then
!          if this is a unit cell, identify this one as a neighbor
!         candidate:
          if(np(node,cell).eq.1)then
            if(separation(cell,i,r_i,1).lt.0)then
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
            if(separation(cell,i,r_cell,0).gt.0)then
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
        include 'sph_common.h'
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

      function separation(cell,i,r,metric)
      integer cell
      include 'sph_common.h'
      d_i_c = 0
c     quadratical metric
      if(metric.eq.1)then
        do k=1,dim
          a = x_cell(node,cell,k) - xp(i,k)
          d_i_c = d_i_c + a*a
        end do
        separation=d_i_c-r*r
      else
c     metric of the maximum:
        do k=1,dim
          a = abs( x_cell(node,cell,k) - xp(i,k) )
          d_i_c = amax(d_i_c,a)
        end do
        separation=d_i_c-r
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
      include 'sph_common.h'
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
            if(separation(cell,i,r,1).lt.0)then
              if(neighb(i).lt.maxneighb)then
                neighb(i) = neighb(i) + 1
                neighb_list(i,neighb(i)) = j
              endif
            endif
          else
            if(separation(cell,i,r_cell,0).gt.0)then
              next = down(node,cell)
              permit(next) = .false.
            endif
          endif
        endif
      end do
      end
