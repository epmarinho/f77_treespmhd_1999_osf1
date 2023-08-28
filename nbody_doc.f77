# 1 "SPH_common.F"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      SYSTEM DEFAULTS FOR THE SPH-TREECODE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      Author:  ERALDO PEREIRA MARINHO.
!
!                  Instituto Astronomico e Geofisico - USP.
!
!                  (011) 577 8599
!
!                  eraldo@astro1.iagusp.usp.br
!
!                  eraldo@iag.usp.ansp.br
!
!                  IAGUSP::ERALDO
!
!      Version: 4.0 10-31-1995
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!      data-type declaration:
!
!
      integer dim,one,two,three,four,eight,down,t,count,row,node,ntbins
     &             ,time_level,global_time_level,timebin,p_list,label
     &             ,deepest,nn,maxoctants,nnode,maxneighb,ichoice,is,ns
     &             ,lst,ifile,ni,vert,ini,nit
     &             ,kwsc,n_timebins,ntbins2,n,next,maxnode,n_p,np,ifirst
     &             ,iiseed,kronecker,levi_civita,sttmnt_count








      real     mass,cellmass,epsilon,eps2_min,eps2_max,g,arest,pi
     &            ,dt_old,p_inertia,x_cell,dt_elapsed,dt,dtm
     &            ,angmom,s,v,xp,theta,theta2,tot_mass
     &            ,E_tot,W_tot,T_tot,U_tot,E_min,E_max,Eo,courant,vmoment,gnode
















      logical  nextnode,permit,convergence
     &            ,you_can,starting
     &            ,converging,ok,nonzero,nonover,nz,nv,ovc
     &            ,not_yet,first



      character esc,cr,cutl*4,rev*4,flsh*4,hltd*4,norm*4,home*4,cls*4
     &              ,answer
      character*132 evolution,name,nameout,namecore,relat,percent







!
!      global parameters:
!
!
      parameter(pi=3.141592654,two_third=.6666667)





      parameter(xnorm=.3183098861)

      parameter(root_3=1.732050808)
      parameter (cc=0.744438)
      parameter (fine=0.56)
!
!            stars parameters:







      parameter(T1=3.37,T2=4.34,ucc=1.343295)
!
!
!      system array sizes:
!      -------------------







      parameter(dim=3,vert=8)









      parameter(nn=32768,maxoctants=262136,nnode=32767)












!
!      usual integers and maximum number of neighbors:
!      -----------------------------------------------
      parameter(zero=0,one=1,two=2,three=3,four=4,eight=8)










!      time resolution:
!      ----------------










      parameter(ntbins=6)

































!      COMMON VARIABLES

!            PART I: STRINGS

!      File names:

      common /answer/ answer
      common /chr/ esc,cr,cutl,rev,flsh,hltd,norm,home,cls
      common /names/       
     &                   evolution
     &                  ,name
     &                  ,nameout
     &                  ,namecore
     &                  ,relat
     &                  ,percent



!
!            PART II: 32-BITS SCALARS




!
      common /choice/                 ichoice
      common /ifnewflag/              first,you_can,not_yet
      common /relate/                 is,ns,ifile,ni,ini,nit,iter
     &                                ,itercounter,iterplus,iter_neighb
      common /nextnode/               nextnode
      common /converge/
     &                                 convergence
     &                               ,converging,ok,nonzero,nonover,nz,nv,ovc
      common /starting_now/             starting
      common /well_sep_node_counter/        kwsc
      common /time_levels/             time_level,global_time_level
      common /time_bins/
     &                               timebin
     &                              ,n_timebins
     &                              ,ntbins2
     &                              ,deepest
     &                              ,dtm












      common /extras2/theta,epsilon,theta2,eps2_max,eps2_min

      common /node/                   node,maxnode
      common /next/                   next
      common /total_number/           n
      common /total_mass/             tot_mass
      common /energies/              E_tot,E_min,E_max,Eo,W_tot
     &                              ,T_tot,U_tot,W_s






      common /Courant_Factor/             courant




!
!
!
!            PART III: SINGLE 3-D VECTORS
!
!
      common /rows/ row (two)
      common /momenta/         vmoment(dim)
     &                        ,angmom (dim)



      common /identity/       kronecker(dim,dim)
      common /levi/           levi_civita(dim,dim,dim)
      common /ifirst/         ifirst      (vert)
      common /count/          count      (vert)
!
!
!            PART IV: Rn-VECTORS
!




      common /particle_mass/ mass(nn)
      common /particle_size/ s(nn)





























      common /tree_seed/             iiseed
!
!       data configuration:
!          nn=n
!            nnode = n - 1
!            maxoctants = vert * nnode
!
!        the array t has the following aspect:
!      -------------------------------------
!                    /t/
!      -------------------------------------
!      /  row  /|           j
!      -------------------------------------
!            i  |  1  |  2  |   3  |   4   
!      -------------------------------------
!            1  |  t  |label| cell | iflg  
!            2  |  "  | "   |  "   |   "   
!            3  |  "  | "   |  "   |   "   
!           ... | ... |...  | ...  |  ...  
!            n  |  "  | "   |  "   |   "   
!
!      where, t = vert * t + cell is a peano parameter, which 
!      distinguish a node from other:
!             cell is the octant number in {1,2,...,vert};
!             iflg = 1 if particle is a leaf, otherwise it equals 0;
!             label is the fixed particle index.
!
      common /t/ t(nn,four)
!
!
!            PART V: Rmn-VECTORS
!










!
!
!            PART VI: Rn3 VECTORS
!
      common /forces/













     &                         g    (nn,dim)

!
!
!            PART VII: Rn nt VECTORS
!
!
      common /time_bin/       n_p(ntbins)
      common /timelist/       p_list(nn,ntbins)
      common /individual_time/dt_old(nn)
      common /time_steps/     dt(0:ntbins)
      common /elapsed_time/   dt_elapsed(ntbins)



!
!            PART VIII: TREE-STRUCTURES  
!
      common /octlist/ lst(nn,vert)
      common /permit/ permit(nnode)
      common /tree_real_0/
     &                 arest        (nnode)
      common /tree_real_1/
     &                 cellmass (nnode, vert)
      common /tree_real_3/
     &                 x_cell   (nnode, vert, dim)
      common /tree_real_6/
     &                 p_inertia(nnode, vert, dim, dim)
      common /tree_intg/
     &                 np       (nnode, vert)
     &                ,label    (nnode, vert)
     &                ,down     (nnode, vert)

      common /gnodes/ gnode(maxoctants,dim)


      common /positions/ xp(nn,dim)
      common /velocities/ v(nn,dim)











      common /boxy_coordinates/box      (nn,dim)

# 1 "SPH_advance.F"
# 104 "SPH_advance.F"



      subroutine Move_particles
      INCLUDE 'SPH_common.h'
      do j=1,dim
        do i_=1,n_p(timebin)
          i=p_list(i_,timebin)
          dt_a =  .5*(dt(timebin)+dt_old(i))
          dt_b = .25*(dt(timebin)-dt_old(i))

          xp(i,j)=xp(i,j)+(v(i,j)+g(i,j)*dt_b)*dt_a



        end do
      end do
      end



      subroutine advance_velocities
      INCLUDE 'SPH_common.h'

      print*,'advance_v'


      call maketree
      call treegravity











      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        do j=1,dim
          v(i,j) = v(i,j) + g(i,j) * dt(timebin)
        end do
      end do














      print*,'advance_v: done'

      end


# 264 "SPH_advance.F"


# 332 "SPH_advance.F"

# 1 "SPH_control.F"




      subroutine Settings
      INCLUDE 'SPH_common.h'
      iiseed=417253423



























      courant = .5
      W_p = 0

      end

      subroutine TimeSet
      INCLUDE 'SPH_common.h'
      dtm=dt(0)
      do j=1,ntbins
        dtm=.5*dtm
      end do

      print*,'dtm =',dtm

      dtm_safe=dtm
      end

      subroutine DoSimulation
      INCLUDE 'SPH_common.h'
      external max, min, amax, amin
      common/general_counter/ icounter
      icounter=0

      print*,'Movie simulation'
      if(global_time_level.eq.1)call Movie (global_time_level)

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









!     perform the simulation now:
      do time_level=1,ninew

          global_time_level = global_time_level + 1

          rtime = global_time_level*dt(0)


          write(*,'(a6,i5,a2,f15.7,a5)')
     &            'Step',time_level,' (',rtime,' [t])'






          if(dt(0).eq.0.0)print*,'simula: warning: dt(0) is null'








          call Integration_scheme (0)
          if(dt(0).eq.0.0)then
            call finalize
            stop
          endif






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




     &            ,W_tot+W_s+2*T_tot     ! 6
     &            ,0.0e+00               ! 7

     &            ,(vmoment(j),j=1,dim)  ! 8..10

     &            ,(angmom(j), j=1,dim)  ! 11..13




          close(2)












!          definitive output:

          call Movie (global_time_level)

          if(is.eq.ns)then
            call finalize
            is=1
            ifile = ifile + 1
          else
            is=is+1
          endif


!         security saving:
!         ------------------
          call writecore


      end do

      print*
      print*,'Simulation done! Have a nice paper!'
      print*

      end

# 1 "SPH_gravity.F"

! 
!                  TREE-GRAVITY SUBROUTINES
! 


      subroutine treepotential
      integer cell,i
      INCLUDE 'SPH_common.h'

      print*,'* treepotential'

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

      print*,'* treepotential end'

      end


      subroutine potential(i,cell)
      integer cell
      INCLUDE 'SPH_common.h'
      common /auxiliar/ aux(nn,dim)
      equivalence (aux, ps)
      dimension ps(nn,dim)
! 
!      Part I. Searching for a well separated cell (good cell)
! 
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
! 
!      Part II. Perform cell contribution to collective potential 
! 



      ep2 = s(label(node,cell))
      ep2 = ep2*ep2 + amax(eps_min,amin(s2,eps2_max))

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

      print*,'* treegravity'

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

      print*,'* treegravity end'

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













        r2    = r2  +  D_x(1,j) * D_x(1,j)
      end do
!     if this is not a good cell, then jump it out; otherwise, set this
!     one as a terminal cell:
      if(np(node,cell).gt.1)then
        !** itheta=-2*int(s(i)/sqrt(r2))
        !** theta_tmp=theta2*(2**itheta)

        !** print*,'theta_tmp=',theta_tmp,' itheta=',itheta

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



      ep2 = s(label(node,cell))
      ep2 = ep2*ep2 + amax(eps_min,amin(s2,eps2_max))

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


# 404 "SPH_gravity.F"

# 1 "SPH_integral.F"
! 
! 
! 
      subroutine Integral_quantities
      INCLUDE 'SPH_common.h'



      call Treepotential
      E_tot = W_tot


      call Angular_momentum

      call Kinetic_energy(1.)
      E_tot = E_tot + T_tot









      end


      subroutine Clear_Ang_Momentum
      INCLUDE 'SPH_common.h'
      do j=1,dim
        angmom(j)=0
      end do
      end

      subroutine Angular_momentum
      INCLUDE 'SPH_common.h'
      call Clear_Ang_Momentum
      do j=1,dim
        vmoment(j) = 0
!        for all system particles do:
        do i=1,n



          v_j = v(i,j)-.5*g(i,j)*dt_old(i)

          vmoment(j) = vmoment(j) + v_j*mass(i)
!          calculate total angular-momentum:
          do k=1,dim
            do l=1,dim
              x_l = xp(i,l)
              angmom(k)=angmom(k) + mass(i)*levi_civita(k,l,j)*x_l*v_j
            end do
          end do
        end do
      end do
      end


# 101 "SPH_integral.F"


      subroutine Kinetic_energy(factor)
      !Don't ask me why I have placed the magnetic energy calculations here!
      INCLUDE 'SPH_common.h'
      ! clear kinetic-energy accumulator:
      T_tot = 0
      do j=1,dim
        do i=1,n

          v_j = v(i,j)-.5*g(i,j)*dt_old(i)



          v2 = v_j * v_j
          T_tot = T_tot + .5*v2*mass(i)
        end do
      end do
      end


















# 1 "SPH_intschm.F"





      subroutine Predictor_corrector (dE)
!      Version Jul  7 1995
!      Revised Sep  5 1995
!              Sep 11 1995
!              Sep 15 1995
!              Jul 11 1996
!            Aug 11 1996
!            Aug 12 1996
!            Aug 13 1996
!            Sep  7 1996
      INCLUDE 'SPH_common.h'




      call move_particles

      call advance_velocities





























      call save_old_timesteps

      print*,'predictor_corrector: done'

      end













# 262 "SPH_intschm.F"



      subroutine Initialize
      INCLUDE 'SPH_common.h'

      print*,'not_yet =',not_yet

      if(not_yet)then

        print*,'Initialize: start'




        n_p(1)      = n
        n_timebins  = 1
        timebin     = 1
        do i=1,n_p(timebin)
          p_list(i,timebin)=i













        end do
        call MakeTree




        you_can = .true.
        not_yet  = .false.
        first    = .false.
        if(dt(0).eq.0.0)return




      endif
      if(ichoice.eq.1)then







        call treegravity




        call Start_list
        do i=1,n
          dt_old(i) = 0
        end do
        ichoice=0
        call predictor_corrector(0)
      endif
      call Start_list




      print*,'initialize: done'

      end



      subroutine Integration_scheme (dE)
      INCLUDE 'SPH_common.h'

      print*,'Integration_scheme'




      call Initialize
      if(dt(0).eq.0.0)return
      call multi_leapfrog (dE)

      call show_bins

      call integral_quantities
      call recover_dtm

      print*,'Integration_scheme: done'

      end
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

      print*
      print*,'writing over ',namecore

      open(3,file=namecore,status='unknown',form='unformatted')

      print*,'successfully opening ',namecore

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
     &             epsilon
      read(1)
     &             E_max
      read(1) 
     &             Eo
      read(1) 
     &             E_min










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


























      close(1)
      theta2=theta*theta
      eps2_max=epsilon*epsilon
      eps2_min=eps2_max*float(n)**(-.3333333)
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
      print *,'enter tree-gravity parameters:'
      print *
      print *,'   aperture parameter:'
      read(*,'(a132)')lntxt
      read(*,*)theta
      if(theta.gt.1)then
        print *,'?: error: invalid theta; must be less than 1.'
        stop
      endif






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




















      print *
      print *,'parameters ok!'
      theta2=theta*theta




      call calc_epsilon


!     find a convenient time-setp:
      if(dt(0).lt.0)then

        dtm = epsilon*sqrt(epsilon/Tot_mass)*float(N)**(1./6)
        dt(0)    = 1024.
        do while(dt(0).gt.dtm)
          dt(0) = .5 * dt(0)
        end do




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


      subroutine calc_epsilon
      logical cont
      INCLUDE 'SPH_common.h'
      print*,'calc_epsilon: start'
      gamma_ = -Tot_Mass*Tot_Mass*float(n)**(-.3333333)
      cont = .true.
      epsilon = 0
      call maketree
      call treepotential
      epsilon  = gamma_/W_tot
      do while(cont)
        epsilon_old = epsilon
        w_old = W_tot
        call treepotential
        epsilon  = 2.*gamma_/(W_tot+w_old)
        cont = abs(epsilon-epsilon_old).gt.(epsilon/32.)
      end do
      eps2_max = epsilon*epsilon
      eps2_min = eps2_max*float(n)**(-.3333333)
      print*,'epsilon_max^2 =',eps2_max,'epsilon_min^2 =',eps2_min
      print*,'calc_epsilon: done'
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
      print *,'base name output sequence        : ',nameout(1:l4)
      print *,'number of next output            : ',ini













      print *,'      pure gravitational simulation!'

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

      character*132 pre1






      character chrd*4
      INCLUDE 'SPH_common.h'
      l1  = index(nameout,']')+1
      l2  = index(nameout(l1:132),'.')-1
      l   = l1+l2
      ext = nameout(l:132)
      call numstring(chrd,mchr,ifile)















      pre1 = nameout(1:l-1) 
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



# 664 "SPH_io.F"

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
        chrd = char(idig+48) 
        itmp = itmp/10
      end do
      end



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
      
      l1  = index(nameout,']')+1
      l2  = index(nameout(l1:132),'.')-1
      l   = l1+l2
      ext = nameout(l:132)
      call numstring(chrd,mchr,leg)
      pre4 = nameout(1:l-1) 
      print*,'saving movie file:',pre4
      n_show=n/npoints
      if(n_show.eq.0)n_show=1
      open(4,file=pre4,status='unknown')
      do i=1,n,n_show
























          write(4,'(3e11.3)')
     &           (xp(i,j),j=1,dim)

      end do
      close(4)
      end

# 1 "SPH_kernel.F"

      function sphKernel(r,hp)
                                                            


      INCLUDE 'SPH_common.h'
      q=r/hp
      if(q.lt.1.)then
          w = 1.5 * (q*q) * (.5*q - 1) + 1
      else
          if(q.lt.2.)then
                w = 2. - q
                w = .25 * (w * w * w)
          else
            w = 0
          endif
      endif





      w = xnorm/(hp*hp*hp) * w

      sphKernel = w
      end



      function D_kernel(r,hp)
                                                            




      INCLUDE 'SPH_common.h'
      q=r/hp
      if(q.lt.1.)then
          dw =  q * (.75*q-1.)
      else
          if(q.lt.2.)then
                dw = 2. - q
                dw = - .25 * (dw * dw)
          else
            dw = 0
          endif
      endif
      overh2 = 1./(hp*hp)






      overh4 = overh2 * overh2
      dw = xnorm * overh4 * dw

      D_kernel = 3. * dw
      end



      function spline(q)
      ! spline kernel (monaghan & lattanzio 1985)
      INCLUDE 'SPH_common.h'
      if(q.lt.1.)then
          w = 1.5 * (q*q) * (.5*q - 1) + 1.
      else
          if(q.lt.2.)then
                w = 2. - q
                w = .25 * (w * w * w)
          else
            w = 0
          endif
      endif
      spline = w
      end


# 1 "SPH_main.F"
      PROGRAM sph_Tree_Gravitation
!
!
!      (individualy varying particle time-steps)
!
!
!
!
!
!
!
!
!
!
!      Author:      Eraldo Pereira Marinho
!
!
!
!
!
!
!
!                  Instituto Astronomico e Geofisico - USP.
!
!                  (011) 577 8599 (branch 235)
!
!                  eraldo@astro1.iagusp.usp.br
!
!                  eraldo@iag.usp.ansp.br
!
!                  FPSP::in%IAGUSP::ERALDO
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!      Version: 24.3 28-5-1994  
!
!
!            Last changings: 
!
!
!
!
!
!      Notes:
!
!            Physical unities are given in
!
!            [l] = 1 pc = 3.086 e+18 cm
!
!            [t] = 1 Myr = 3.156 e+13 s
!
!            [m] = 2.224 e+2 solar masses = 442.35 e+33 g
!            
!
!
!
!
! 
!
!23456789012345678901234567890123456789012345678901234567890123456789012
!        1         2         3         4         5         6         7
!
      print*,'////////////////////////////////////////////'
      print*



      print*
      print*,'////////////////////////////////////////////'
      print*
      print*,'Vector/Parallel SPH treecode v. 11-12-1995'
      print*,'      Version 5.0'
      print*
      print*
      print*,'Designed by ERALDO PEREIRA MARINHO'
      print*
      print*,'      Instituto Astronomico e Geofisico.'
      print*,'      Universidade de Sao Paulo - Brazil'
      print*
      print*,'(P) 1992. All Rights Reserved.'
      print*
      print*,'////////////////////////////////////////////'
      print*
      print*
      call dosimulation
      print *
      print *,'Done!'
!
      end
! 

# 1 "SPH_sseeds.F"
# 45 "SPH_sseeds.F"


# 1 "SPH_timing.F"



      subroutine show_bins
      integer iti,cksum
      INCLUDE 'SPH_common.h'
      cksum=0
      print*
      write(*,'(a9)')
     &      'covered '
      do iti=1,n_timebins
        write(*,'(f10.2,x,a10,x,i2,x,i6,x,a10,x,e16.7)')
     &        100*amin(dt_elapsed(iti)/dt(iti-1),1.)
     &        ,'% in level',iti,n_p(iti),' particles',dt(iti)
        cksum=cksum+n_p(iti)
      end do
      print*
      print*,'check-sum:',cksum
      end


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



          !print*,'dt_elapsed=',dt_elapsed(timebin),' timebin=',timebin
      goto 1111
 1114 continue!--> end while
      timebin = timebin - 1
      if(timebin.eq.0)return
      goto 1112!--> recursive return:
      end


      subroutine recover_dtm







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



        gx=abs(g(i,j))

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
















      dt_Courant = ttc

      end


      subroutine save_old_timesteps
      integer il,i
      INCLUDE 'SPH_common.h'
      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        dt_old(i) = dt(timebin)
      end do
      end































# 1 "SPH_tree.F"
! 
!
!      TREE CONSTRUCTION AND TREE-ACCESS SUBROUTINES:
!
! 
!
      subroutine SetBox
      INCLUDE 'SPH_common.h'
      real aux,vertix



      real scale,s1bit,ar,c,a

      common/aux/ aux(nn,dim)
      equivalence (aux,vertix)
      dimension vertix(1,dim)
!
!
!      * take the lower vertix from the trihedron enclosing the system:
!          vertix(j)=min(xp(i,j))
!      * set the major cube arest:
!          arest = max(xp(i,j))-min(xp(i,j))
!
      do j=1,dim
          vertix(1,j)=xp(1,j)
      end do
      c=xp(1,1)
      a=c
      do j=1,dim
          do i=2,n
            if(xp(i,j).gt.c)then
                c=xp(i,j)
            else
                if(xp(i,j).lt.vertix(1,j))vertix(1,j)=xp(i,j)
            endif
          end do
          if(vertix(1,j).lt.a)a=vertix(1,j)
      end do
      ar=c-a
!
!      * determina o cubo circunscrito cujas arestas sao potencias
!        inteiras de 2:
!
      s1bit=1.
      if(ar.gt.1.)then
          do while(s1bit.lt.ar)
            s1bit=s1bit+s1bit
          end do
      else
          ar=2*ar
          do while(s1bit.gt.ar)
            s1bit=0.5*s1bit
          end do
      endif
      if(s1bit.ne.0)then
          scale=1d0/s1bit
          arest(1)=0.5*s1bit
      else
          print *,
     &      '%err: particles were distributed into a singularity!!'
          stop
      endif
!
!      renormaliza as coordenadas das particulas para um cubo
!      unitario:
!
      do j=1,dim
          do i=1,n




            box(i,j)=scale*(xp(i,j)-vertix(1,j))

          end do
      end do
      end
!
!
      subroutine MakeTree
!      this procedure generates a tree data-structure for a given 
!      set of n mass-points in 3-dim space. the method consists on 
!      updating the particles in all corresponding nodes in each level
!      at once, which implies in a one descenting construction, 
!      discarding recursive procedures. the philosophy of tree-
!      descendent construction has the principal advantage on its easy
!      adaptation for vector machines.
!
!    data configuration:
!      nn=n
!            nnode = n - 1
!            maxoctants = eight * nnode
!
!
!
!
      integer i,i1,nsave
      real volume,dens_mean
      INCLUDE 'SPH_common.h'
      common /indices/ i,i1
!
!
!      inicializacao
!

      print*,'maketree: start'


      call SetBox
      nsave=n
      do i=1,n
          t(i,1)=0
          t(i,four)=0
          t(i,2)=i
      end do









      node = 1
      next = node + 1
!
!      while any particle is non-leaf yet, continue building tree:
      do while(n.gt.1)
          call FindOctants(1,n)
          i=1
          do while(i.le.n)
            !make node by sharing cubic-cells
            !with respective particles:
            call MakeNode






!            link node to the tree:
            call LinkTree
          end do
          call NoLeaves
      end do
!
!      finalization of tree-construction:
      maxnode=node-1
      n=nsave

      print*,'Tree completed! n=',n

!
# 187 "SPH_tree.F"


      call Quadrupole


      print*,'maketree: done'

      end 



      subroutine MakeNode
      logical changing
      integer cell,ip,i1,i,j
      INCLUDE 'SPH_common.h'
      common /indices/ i,i1
      common /non_singularity/ isingularity
      isingularity=0
      call CLR
      i1=i
      do while(.not.nextnode)
          ip=t(i,2)
          cell=t(i,3)+1
          count(cell)=count(cell)+1
          cellmass(node,cell)=cellmass(node,cell)+mass(ip)
          lst(count(cell),cell)=ip
!         calculate the first-order position-momentum:
          do j=1,dim
             x_cell(node,cell,j)=x_cell(node,cell,j)+xp(ip,j)*mass(ip)
          end do
          if(count(cell).eq.1)ifirst(cell)=i
          i=i+1
          if(i.le.n)changing=t(i,1).gt.t(i-1,1)
          nextnode=changing.or.(i.gt.n)
          if(nextnode)call NonDegeneracy
      end do
      end 




      subroutine LinkTree
      integer cell,j,l,ip
      INCLUDE 'SPH_common.h'

      do cell=1,vert



        np(node,cell)=count(cell)
!        make centroid:
        if(count(cell).ge.1)then
!          resume cell centroide:
          do j=1,dim
           x_cell(node,cell,j)=x_cell(node,cell,j)/cellmass(node,cell)
          end do
        endif
!
        if(count(cell).gt.1)then
            down(node,cell)=next
            arest(next)=0.5*arest(node)
            next=next+1
        else
            down(node,cell)=0
            if(count(cell).eq.1)then
                l=ifirst(cell)
                ip=t(l,2)
                s(ip)=arest(node)
                t(l,four) = 1
                label(node,cell)=ip
            endif
        endif
      end do
      node=node+1
      if(node.gt.nnode)then
        print*,'LinkTree: ERROR: tree overflow'
        print*,'node=',node,'  nnode=',nnode
        stop
      endif
      end




      subroutine NonDegeneracy
      integer cell,i,i1,idegen
      INCLUDE 'SPH_common.h'
      common /non_singularity/ isingularity
      common /indices/ i,i1
      idegen=0
      do cell=1,vert
          if(count(cell).eq.0)idegen=idegen+1
      end do
      if(idegen.eq.(vert-1))then
          call CLR
          isingularity=isingularity+1
          call FindOctants(i1,i-1)
          arest(node)=0.5*arest(node)
          i=i1
      else
          isingularity=0
      endif
      end 






      subroutine FindOctants(ia,ib)
      integer cell,i,ia,ib,j,ip



      INCLUDE 'SPH_common.h'
      common /non_singularity/ isingularity
      do i=ia,ib
          ip=t(i,2)
          cell=0
          do j=1,dim
                cell=2*cell
            box(ip,j)=2*box(ip,j)


            if(isingularity.gt.23)then
              isingularity=0

              print*,'Tree-singularity occurring: assuming dynamic round-off'
              print*,'old box: ',box(ip,j)




              xperr = -1e-6*(1+ran(iiseed))

              box(ip,j) = box(ip,j)*(1+xperr)

              print*,'new box: ',box(ip,j)
              print*,'x-error: ',xperr

            endif


            if(box(ip,j).gt.1)then
                box(ip,j)=box(ip,j)-1
                cell=cell+1
            endif

            if(box(ip,j).gt.1)then
              print*,'FindOctants: error'
              print*,'box=',box(ip,j)
              stop
            endif

          end do
          t(i,3)=cell
      end do
      end 




      subroutine CLR
      integer cell,j
      INCLUDE 'SPH_common.h'

      do cell=1,vert



          count(cell) = 0
          cellmass(node,cell) = 0
          do j=1,dim
            x_cell(node,cell,j) = 0
          end do
      end do
      nextnode = .false.
      end
!
!
      subroutine NoLeaves
      integer k,i,j
      logical noleaf
!
      INCLUDE 'SPH_common.h'
!
!      copia somente as linhas nao-folha
!      e dispensa a coluna dos "leaf-flags"(4):
!
      k=0
      do i=1,n
          noleaf=(t(i,four).ne.1)
          if(noleaf)then
            k=k+1
            do j=1,four
                t(k,j)=t(i,j)
            end do
          endif
      end do
      n=k
      if(n.eq.0)return
!
!      calcula o novo parametro de peano promovendo os octantes para o 
!      status de no':
!
      do i=1,n

          t(i,1)=vert*t(i,1)+t(i,3)



      end do
!
!      ordena o array t por heap-sort dispensando a coluna dos 
!      octantes (col. 3):
!
      call hsort
!
!      reduz o parametro de peano a fim de evitar estouro de inteiro:
!
      k=0
      do i=2,n
          if(t(i,1).ne.t(i-1,1))k=k+1
          t(i-1,3)=k
      end do
      do i=2,n
          t(i,1)=t(i-1,3)
      end do
      t(1,1)=0
!
      end 
!
!
      subroutine hsort
!
!            method used: heapsort (num. rec. 1986)
!
      integer l,ir,k,i,j
      INCLUDE 'SPH_common.h'
!
!
      l=n/2+1
      ir=n
   10      continue
          if(l.gt.1)then
            l=l-1
            do k=1,2
                row(k)=t(l,k)
                t(l,k)=t(1,k)
            end do
          else
            do k=1,2
                row(k)=t(ir,k)
                t(ir,k)=t(1,k)
            end do
            ir=ir-1
            if(ir.eq.1)then
                do k=1,2
                  t(1,k)=row(k)
                end do
                return
            endif
          endif
          i=l
          j=l+l
   20          if(j.le.ir)then
            if(j.lt.ir)then
                if(t(j,1).lt.t(j+1,1))j=j+1
            endif
            if(row(1).lt.t(j,1))then
                do k=1,2
                  t(i,k)=t(j,k)
                end do
                i=j
                j=j+j
            else
                j=ir+1
            endif
          go to 20
          endif
          do k=1,2
            t(i,k)=row(k)
          end do
      go to 10
      end
!
!
      subroutine forbid
!
      integer cell
!
      INCLUDE 'SPH_common.h'
!



      do cell=1,vert

          next=down(node,cell)
          if(next.gt.0)permit(next)=.false.
      end do
      end
!
!

      subroutine quadrupole
      integer node1,m1,i,j,node2,m2
      real aux,D_x
      INCLUDE 'SPH_common.h'
      common /auxiliar/ aux(nn,dim)
      dimension D_x(nn,dim)
      equivalence (aux, D_x)
!
      do node1=maxnode,1,-1
        do m1=1,eight
!
          do i=1,dim
            do j=i,dim
            p_inertia(node1,m1,i,j)=0
            p_inertia(node1,m1,j,i)=0
            end do
          end do
!
          if(np(node1,m1).gt.1)then
!
            node2  =  down(node1,m1)
!
            do m2=1,eight
!
            do i=1,dim
              D_x(1,i)  =  x_cell(node1,m1,i) - x_cell(node2,m2,i)
            end do
!
!            parallel axes theorem:
            do j=1,dim
                    do i=j,dim
!
                p_inertia(node1,m1,i,j) = p_inertia(node1,m1,i,j) +
     &                                          p_inertia(node2,m2,i,j) +
     &                         D_x(1,i)*cellmass(node2,m2)*D_x(1,j)
!
                p_inertia(node1,m1,j,i) = p_inertia(node1,m1,i,j)
!
                    end do
            end do
!
            end do
!
          endif
!
        end do
      end do
      end





























# 1 "SPH_usr.F"
# 248 "SPH_usr.F"

# 1 "SPH_usr_accel.F"
# 106 "SPH_usr_accel.F"

