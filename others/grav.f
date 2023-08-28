
      subroutine treegravity
      integer cell,il,i,j,l
      INCLUDE 'SPH_common.h'
      !dimension xa(dim),ga(dim)
      print*,'treegravity'
      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        do node=1,maxnode
          permit(node)=.true. 
        end do
        kwsc=0
        ! Performing tree-walk
        do node=1,maxnode
          ! vectorization is broken here:
          if(permit(node))then
            do cell=1,eight
              if(np(node,cell).gt.0)call gravity(i,cell)
            end do
          else
            call forbid! cell descent
          endif
        end do
        !-----------------------------
        !do j=1,dim
        !  xa(j)=xp(i,j)
        !end do
        !call phi_galax(xa,ga,epsilon)
        ! Now, the particle has a list
        ! of kwsc cell'contributions:
        !-----------------------------
        do j=1,dim
          g(i,j) = 0
          do l=1,kwsc
            g(i,j)=g(i,j)+gnode(l,j)
            !g(i,j)=g(i,j)+ga(j)
          end do
        end do
        !-----------------------------
      end do
      print*,'done!'
      return
      end


      subroutine gravity(i,cell)
      integer cell
      INCLUDE 'SPH_common.h'
      dimension D_x(dim)
!     square cell size (it's easy with a single bit value):
      s2 = arest(node)*arest(node)
!     calculate cell-to-particle squared hard-distance.
!     r2 <-- r * r
      r2=0
      do j=1,dim
        ! D_x is the position vector of the particle with respect to
        ! the observed cell:
        !
        D_x(j) = xp(i,j) - x_cell(node,cell,j)
        r2    = r2  +  D_x(j) * D_x(j)
      end do
!     if this is not a good cell, then jump it out; otherwise, set this
!     one as a terminal cell:
      if(np(node,cell).gt.1)then
        if((theta2*r2).lt.s2)then
          return
        else
          permit(down(node,cell))=.false.
        endif
      else
        if(label(node,cell).eq.i)return
        !since we may not calculate particle's self gravity:
      endif
!     count the good-cell entry:
      kwsc = kwsc + 1
!     calculate softened cell-to-particle linear, cubic and squared
!     distances respectively:
      ep2 = epsilon2
      r2 = r2 + ep2
      r  = sqrt(r2)
      r3 = r2 * r
      do j=1,dim
        D_x(j)  =  D_x(j) / r2
      end do
      phi0 = - cellmass(node,cell) / r
      do j=1,dim
        gnode(kwsc,j) = phi0 * D_x(j)
      end do
      diag=0
      trace=0
      do k=1,dim
          diag = diag + D_x(k)*p_inertia(node,cell,k,k)*D_x(k)
          trace = trace + p_inertia(node,cell,k,k)
      end do
      triang=0
      do k=2,dim
          do l=1,k-1
            triang=triang + D_x(k)*p_inertia(node,cell,k,l)*D_x(l)
          end do
      end do
      doub = diag + 2 * triang
      phi2 = -.5*(3*doub - trace/r2)/r
      phi   =   phi0   +    phi2
      do j=1,dim
        singl=0
        do k=1,dim
            singl = singl + p_inertia(node,cell,j,k)*D_x(k)
        end do
        gnode(kwsc,j) = gnode(kwsc,j) + 
     &        5 * phi2 * D_x(j) + (  3*singl - D_x(j) * trace  ) / r3
      end do
      end
