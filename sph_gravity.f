#ifndef _NON_SELFGRAVITATING
!///////////////////////////////////////////////////////////////////////
!                  TREE-GRAVITY SUBROUTINES
!///////////////////////////////////////////////////////////////////////

#  ifndef _TUBE
      subroutine treepotential
      integer cell,i
      INCLUDE 'SPH_common.h'
#ifdef _VERBOSE
      print*,'* treepotential'
#endif
      W_tot = 0
      W_s   = 0
      do i=1,n
        do node=1,maxnode
          permit(node)=.true.
        end do
        do node=1,maxnode
          !this descent-control breaks vectorization
          if(permit(node))then
            do cell=1,eight
              if(np(node,cell).gt.0)call potential(i,cell)
            end do
          else
            call forbid
          endif
        end do
      end do
#ifdef _VERBOSE
      print*,'* treepotential end'
#endif
      end


      subroutine potential(i,cell)
      integer cell
      INCLUDE 'SPH_common.h'
      common /auxiliar/ aux(nn,dim)
      equivalence (aux, ps)
      dimension ps(nn,dim)
!///////////////////////////////////////////////////////////////////////
!      Part I. Searching for a well separated cell (good cell)
!///////////////////////////////////////////////////////////////////////
!      square up the box size:
      s2 = arest(node)
      s2=s2*s2
!      calculate cell-to-particle actual squared distance.
!      r2 <-- r * r
      r2=0
      do j=1,dim
        ps(1,j) = xp(i,j) - x_cell(node,cell,j)
        r2    = r2  +  ps(1,j) * ps(1,j)
      end do
!      if it is not a good cell then go out
      if(np(node,cell).gt.1)then
        if((theta2*r2).lt.s2)then
          return
        else
          permit(down(node,cell))=.false.
        endif
      else
        if(label(node,cell).eq.i)return
      endif
!///////////////////////////////////////////////////////////////////////
!      Part II. Perform cell contribution to collective potential 
!///////////////////////////////////////////////////////////////////////
#    if defined (_CONSTANT_SOFTENING)
      ep2=eps2_max
#    else
      ep2 = s(label(node,cell))
      ep2 = ep2*ep2 + amax(eps2_min,amin(s2,eps2_max))
#    endif
!     square softened distance:
      r2 = r2 + ep2
!     softened distance:
      r  = sqrt(r2)
!     cubic softened distance:
      r3 = r2 * r
!     calculate pole-vector:
      do j=1,dim
        ps(1,j)=ps(1,j)/r2
      end do
!     --------------------------------
!     perform monopole approximation
!     ------------------------------
      phi0 = - cellmass(node,cell) / r
!     --------------------------------
!     perform quadrupole approximation
!     --------------------------------
!     diagonal sumation:
      diag=0
      trace=0
      do k=1,dim
          diag = diag + ps(1,k)*p_inertia(node,cell,k,k)*ps(1,k)
          trace = trace + p_inertia(node,cell,k,k)
      end do
!      triangular sumation: quadrupole is symmetric!
      triang=0
      do k=2,dim
        do l=1,k-1
          triang = triang + ps(1,k)*p_inertia(node,cell,k,l)*ps(1,l)
        end do
      end do
!      doub gets the double internal-product:    ^p . ^Q^ . ^p
      doub = diag + 2 * triang
!      perform quadrupole contribution:
      phi2 = -.5*(3*doub - trace/r2)/r
!      calculate potential contribution
      phi = phi0 + phi2
      E_pair  =  mass(i) * phi
!      add "E_pair" to total potential-energy "W_tot"
      W_tot    =    W_tot    +    .5  *  E_pair
!      calculate soft-virial contribution
      E_s  =  .5  *  mass(i) * cellmass(node,cell) * ep2 / r3
!      add "E_s" to total soft-virial "W_s"
      W_s    =    W_s    +    E_s
      return
      end


      subroutine treegravity
      integer cell,il,i,j,l
      INCLUDE 'SPH_common.h'
      !dimension xa(dim),ga(dim)
#ifdef _VERBOSE
      print*,'* treegravity'
#endif
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
        do j=1,dim
            g(i,j) = 0
            do l=1,kwsc
            g(i,j)=g(i,j)+gnode(l,j)
            !g(i,j)=g(i,j)+ga(j)
            end do
        end do
        !-----------------------------
      end do
#ifdef _VERBOSE
      print*,'* treegravity end'
#endif
      return
      end


      subroutine gravity(i,cell)
      integer cell
      INCLUDE 'SPH_common.h'
      common /auxiliar/ aux(nn,dim)
      dimension D_x(nn,dim)
      equivalence (aux, D_x)
**************************************************************
*     performs a tree-descent by searching for well-separated
*     cells in the tree
**************************************************************
!     square cell size (it's easy with a single bit value):
      s2 = arest(node)*arest(node)
!     calculate cell-to-particle squared hard-distance.
!     r2 <-- r * r
      r2=0
      do j=1,dim
        ! D_x is the position vector of the particle with respect to
        ! the observed cell:
        !
        D_x(1,j) = xp(i,j) - x_cell(node,cell,j)
#    ifdef _SYNCHRONIZE_
        !*************************************************************
        ! inserted on Jul 17 1995, but it did not solved the problem
        ! of numerical cooling:
        !*************************************************************
        if(np(node,cell).eq.1)then
          l=label(node,cell)
          ddt=dt_elapsed(timebin)-dt_elapsed(p_home(l))
          dx=ddt*(v(l,j)-.5*ddt*accel(l,j))
          D_x(1,j) = D_x(1,j) - dx
        endif
        !*************************************************************
#    endif
        r2    = r2  +  D_x(1,j) * D_x(1,j)
      end do
!     if this is not a good cell, then jump it out; otherwise, set this
!     one as a terminal cell:
      if(np(node,cell).gt.1)then
        !** itheta=-2*int(s(i)/sqrt(r2))
        !** theta_tmp=theta2*(2**itheta)
#    ifdef _VERBOSE
        !** print*,'theta_tmp=',theta_tmp,' itheta=',itheta
#    endif _VERBOSE
        !** if((theta_tmp*r2).lt.s2)then
        if((theta2*r2).lt.s2)then
          return
        else
          permit(down(node,cell))=.false.
        endif
      else
        if(label(node,cell).eq.i)return
        !since we may not calculate particle's self gravity:
      endif
**************************************************************
!     count the good-cell entry:
      kwsc = kwsc + 1
!     calculate softened cell-to-particle linear, cubic and squared
!     distances respectively:
#    if defined (_CONSTANT_SOFTENING)
      ep2=eps2_max
#    else
      ep2 = s(label(node,cell))
      ep2 = ep2*ep2 + amax(eps2_min,amin(s2,eps2_max))
#    endif
      r2 = r2 + ep2
      r  = sqrt(r2)
      r3 = r2 * r
!      calculate "pole-vector" (softened unit vector):
      do j=1,dim
        D_x(1,j)  =  D_x(1,j) / r2
      end do

! **************************
! perform cell contributions
! **************************

!     ******************************
!     perform monopole approximation
!     ******************************

!     potential:
      phi0 = - cellmass(node,cell) / r
      !for testing:
      !*** if(phi0.eq.0)then
            !*** !print*,'gravity: error: cellmass is zero'
            !*** stop
      !*** endif

!      monopole-gravity-field:
      do j=1,dim
        gnode(kwsc,j) = phi0 * D_x(1,j)
      end do

!     ********************************
!     perform quadrupole approximation
!     ********************************

!     diagonal sumation:
      diag=0
      trace=0
      do k=1,dim
          diag = diag + D_x(1,k)*p_inertia(node,cell,k,k)*D_x(1,k)
          trace = trace + p_inertia(node,cell,k,k)
      end do
!      triangular sumation: quadrupole is symetric!
      triang=0
      do k=2,dim
          do l=1,k-1
            triang=triang + D_x(1,k)*p_inertia(node,cell,k,l)*D_x(1,l)
          end do
      end do
!      doub gets double the internal-product:    ^p . ^q^ . ^p
      doub = diag + 2 * triang
!      perform quadrupole cell contribution to potential:
      phi2 = -.5*(3*doub - trace/r2)/r
!      calculate potential contribution
!      total = monopole + quadrupole
      phi   =   phi0   +    phi2

      do j=1,dim
!        singl gets the single internal-product:    ^q^ . ^p
        singl=0
        do k=1,dim
            singl = singl + p_inertia(node,cell,j,k)*D_x(1,k)
        end do
        gnode(kwsc,j) = gnode(kwsc,j) + 
     &        5 * phi2 * D_x(1,j) + (  3*singl - D_x(1,j) * trace  ) / r3
      end do
      return
      end
#  endif

#elif defined(_SIMPLER_DISK_MODEL)

#  include "corepot.f"

#elif defined(_GALAXY_2D)
      subroutine extern_gravity_2d
      INCLUDE 'SPH_common.h'
      !!print*,'extern_gravity_2d'

      ! physical units must be given as [l]=1 kpc, [t]=100 Myr
      ! unperturbed parameters:
      v1=33.77 != 330 km/s
      v2=86.98 != 850 km/s
      ra1=18.0 !kpc
      rb1=1.50 !kpc
      ra2=0.73 !kpc
      rb2=0.30 !kpc
      rcore=0.1!kpc
      rcore=rcore*rcore

      ! perturbation parameters:
      pitch=-(7.0/180.0)*pi
      pitch=1/tan(pitch)
      R0=7.5 !kpc
      Ampl=-2.612 ! [v^2]
      r=R0
      Omega_p=2.4291!(v1*exp(-r/ra1-rb1/r)+v2*exp(-r/ra2-rb2/r))/r
      !!print*,'Omega_p=',Omega_p
      nspirals=2

      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        !!print*,'i=',i
        r2=0
        do j=1,dim
          r2=r2+xp(i,j)*xp(i,j)
        end do
        r=sqrt(r2+rcore)
        !!print*,'r=',r

        !*** unperturbed potential ***
        vc=v1*exp(-r/ra1-rb1/r)+v2*exp(-r/ra2-rb2/r)
        !!print*,'vc=',vc
        gc=-vc*vc/r
        !!print*,'g=',gc

#  ifdef _SPIRAL
        !*** spiral perturbation ***
        angle=acos(xp(i,1)/r)
        if(xp(i,2).lt.0)angle=-angle
        gs=Ampl
#    if defined(_SLOW_SWITCH_ON)
     &    *(1-exp(-.2*Omega_p*time(timebin)))
#    endif
     &    *sin(
     &      nspirals*(
     &        pitch*log(r/R0)-angle
#    if !defined(_ROTATING_FRAME)
     &        +Omega_p*time(timebin)
#    endif
     &      )
     &    )
     &    /r
        g(i,1)=gs*(pitch*xp(i,1)+xp(i,2))/r
        g(i,2)=gs*(pitch*xp(i,2)-xp(i,1))/r
#  endif _SPIRAL

        !*** finalization
        do j=1,dim
#  ifdef _SPIRAL
          g(i,j)=g(i,j)+gc*xp(i,j)/r
#  else
          g(i,j)=gc*xp(i,j)/r
#  endif _SPIRAL
#  ifdef _ROTATING_FRAME
          !*** include centrifugal force (don't forget the Coriolis one!)
          g(i,j)=g(i,j)+Omega_p**2*xp(i,j)
#  endif _ROTATING_FRAME
        end do

      end do

#  ifdef _VERBOSE
      write(*,'(a16,e15.7,a11,i2)')
     &'potential time=',time(timebin),' time-bin=',timebin
#  endif
      end

#  ifdef _ROTATING_FRAME
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
#  endif _ROTATING_FRAME

#endif _NON_SELFGRAVITATING
