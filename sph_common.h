# 1 "SPH_common.F"
# 1 "<built-in>"
# 1 "<command line>"
# 1 "SPH_common.F"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SYSTEM DEFAULTS FOR THE SPH-TREECODE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Author: ERALDO PEREIRA MARINHO.
!
! Instituto Astronomico e Geofisico - USP.
!
! (011) 577 8599
!
! eraldo@astro1.iagusp.usp.br
!
! eraldo@iag.usp.ansp.br
!
! IAGUSP::ERALDO
!
! Version: 4.0 10-31-1995
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
! data-type declaration:
!
!
      integer dim,one,two,three,four,eight,down,t,count,row,node,ntbins
     & ,time_level,global_time_level,timebin,p_list,label
     & ,deepest,maxoctants,nnode,maxneighb,ichoice,is,ns
     & ,lst,ifile,ni,vert,ini,nit
     & ,kwsc,n_timebins,ntbins2,n,next,maxnode,n_p,np,ifirst
     & ,iiseed,kronecker,levi_civita,sttmnt_count
# 40 "SPH_common.F"
      real mass,cellmass,epsilon,eps2_min,eps2_max,g,arest,pi
     & ,dt_old,p_inertia,x_cell,dt_elapsed,dt,dtm
     & ,angmom,s,v,xp,theta,theta2,tot_mass
     & ,E_tot,W_tot,T_tot,U_tot,E_min,E_max,Eo,courant,vmoment,gnode
# 60 "SPH_common.F"
      logical nextnode,permit,convergence
     & ,you_can,starting
     & ,converging,ok,nonzero,nonover,nz,nv,ovc
     & ,not_yet,first



      character esc,cr,cutl*4,rev*4,flsh*4,hltd*4,norm*4,home*4,cls*4
     & ,answer
      character*132 evolution,name,nameout,namecore,relat,percent







!
! global parameters:
!
!
      parameter(pi=3.141592654,two_third=.6666667)





      parameter(xnorm=.3183098861)

      parameter(root_3=1.732050808)
      parameter (cc=0.744438)
      parameter (fine=0.56)
!
! stars parameters:







      parameter(T1=3.37,T2=4.34,ucc=1.343295)
!
!
! system array sizes:
! -------------------







      parameter(dim=3,vert=8)
# 133 "SPH_common.F"
      parameter(nn=262144,maxoctants=2097144,nnode=262143)


!
! usual integers and maximum number of neighbors:
! -----------------------------------------------
      parameter(zero=0,one=1,two=2,three=3,four=4,eight=8)
# 150 "SPH_common.F"
! time resolution:
! ----------------
# 162 "SPH_common.F"
      parameter(ntbins=6)
# 196 "SPH_common.F"
! COMMON VARIABLES

! PART I: STRINGS

! File names:

      common /answer/ answer
      common /chr/ esc,cr,cutl,rev,flsh,hltd,norm,home,cls
      common /names/
     & evolution
     & ,name
     & ,nameout
     & ,namecore
     & ,relat
     & ,percent



!
! PART II: 32-BITS SCALARS




!
      common /choice/ ichoice
      common /ifnewflag/ first,you_can,not_yet
      common /relate/ is,ns,ifile,ni,ini,nit,iter
     & ,itercounter,iterplus,iter_neighb
      common /nextnode/ nextnode
      common /converge/
     & convergence
     & ,converging,ok,nonzero,nonover,nz,nv,ovc
      common /starting_now/ starting
      common /well_sep_node_counter/ kwsc
      common /time_levels/ time_level,global_time_level
      common /time_bins/
     & timebin
     & ,n_timebins
     & ,ntbins2
     & ,deepest
     & ,dtm
# 250 "SPH_common.F"
      common /extras2/theta,epsilon,theta2,eps2_max,eps2_min

      common /node/ node,maxnode
      common /next/ next
      common /total_number/ n
      common /total_mass/ tot_mass
      common /energies/ E_tot,E_min,E_max,Eo,W_tot
     & ,T_tot,U_tot,W_s






      common /Courant_Factor/ courant




!
!
!
! PART III: SINGLE 3-D VECTORS
!
!
      common /rows/ row (two)
      common /momenta/ vmoment(dim)
     & ,angmom (dim)



      common /identity/ kronecker(dim,dim)
      common /levi/ levi_civita(dim,dim,dim)
      common /ifirst/ ifirst (vert)
      common /count/ count (vert)
!
!
! PART IV: Rn-VECTORS
!




      common /particle_mass/ mass(nn)
      common /particle_size/ s(nn)
# 324 "SPH_common.F"
      common /tree_seed/ iiseed
!
! data configuration:
! nn=n
! nnode = n - 1
! maxoctants = vert * nnode
!
! the array t has the following aspect:
! -------------------------------------
! /t/
! -------------------------------------
! / row /| j
! -------------------------------------
! i | 1 | 2 | 3 | 4
! -------------------------------------
! 1 | t |label| cell | iflg
! 2 | "  | " | "   |   "
! 3 | "  | " | "   |   "
! ... | ... |... | ... | ...
! n | "  | " | "   |   "
!
! where, t = vert * t + cell is a peano parameter, which
! distinguish a node from other:
! cell is the octant number in {1,2,...,vert};
! iflg = 1 if particle is a leaf, otherwise it equals 0;
! label is the fixed particle index.
!
      common /t/ t(nn,four)
!
!
! PART V: Rmn-VECTORS
!
# 366 "SPH_common.F"
!
!
! PART VI: Rn3 VECTORS
!
      common /forces/
# 384 "SPH_common.F"
     & g (nn,dim)

!
!
! PART VII: Rn nt VECTORS
!
!
      common /time_bin/ n_p(ntbins)
      common /timelist/ p_list(nn,ntbins)
      common /individual_time/dt_old(nn)
      common /time_steps/ dt(0:ntbins)
      common /elapsed_time/ dt_elapsed(ntbins)



!
! PART VIII: TREE-STRUCTURES
!
      common /octlist/ lst(nn,vert)
      common /permit/ permit(nnode)
      common /tree_real_0/
     & arest (nnode)
      common /tree_real_1/
     & cellmass (nnode, vert)
      common /tree_real_3/
     & x_cell (nnode, vert, dim)
      common /tree_real_6/
     & p_inertia(nnode, vert, dim, dim)
      common /tree_intg/
     & np (nnode, vert)
     & ,label (nnode, vert)
     & ,down (nnode, vert)

      common /gnodes/ gnode(maxoctants,dim)


      common /positions/ xp(nn,dim)
      common /velocities/ v(nn,dim)
# 433 "SPH_common.F"
      common /boxy_coordinates/box (nn,dim)
