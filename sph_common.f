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
     &             ,deepest,maxoctants,nnode,maxneighb,ichoice,is,ns
     &             ,lst,ifile,ni,vert,ini,nit
     &             ,kwsc,n_timebins,ntbins2,n,next,maxnode,n_p,np,ifirst
     &             ,iiseed,kronecker,levi_civita,sttmnt_count
#ifndef _N_BODY
     &             ,iter,itercounter,iterplus,iter_neighb
     &             ,n_fix,neighb_list,neighb
     &             ,iseed,jseed,kseed,lseed,mseed,nseed,tol_nfix
#endif
#ifdef _TUBE
      real *8 box
#endif
      real     mass,cellmass,epsilon,eps2_min,eps2_max,g,arest,pi
     &            ,dt_old,p_inertia,x_cell,dt_elapsed,dt,dtm
     &            ,angmom,s,v,xp,theta,theta2,tot_mass
     &            ,E_tot,W_tot,T_tot,U_tot,E_min,E_max,Eo,courant,vmoment,gnode
#ifndef _N_BODY
     &            ,lambda,luminosity,two_third,xnorm,root_3,sc,W_p,W_s
     &            ,cc,accel,rho,v_o,v_old
     &            ,u,u_dot,u_o,h,h_dot
     &            ,T1,T2,ucc,u_mass,u_leng
     &            ,fine,u_time,u_dens,u_den_n,XX,YY,ZZ,dtm_safe
     &            ,alpha,beta,eta,guessing
     &            ,P,QQ,P_acc,q_acc,grad_W,u_dot_visc,Temp,cs,u_old,vmonag
#ifndef _SPH_ADIABAT
     &            ,u_lum,u_cool,u_ph,E_n
#if defined(_SPH_STR) || defined(_SPH_STR_SN)
#endif
     &            ,r_sol,SN_lif_time,SN_ene,t_SN,star_mass,n_str,y,str_growth
     &            ,star_min,star_max,SN_mass,solar,u_sol
#endif
#endif
      logical  nextnode,permit,convergence
     &            ,you_can,starting
     &            ,converging,ok,nonzero,nonover,nz,nv,ovc
     &            ,not_yet,first
#ifndef _N_BODY
     &            ,h_list,i_list,gas
#endif
      character esc,cr,cutl*4,rev*4,flsh*4,hltd*4,norm*4,home*4,cls*4
     &              ,answer
      character*132 evolution,name,nameout,namecore,relat,percent
#if defined(_SPH_STR) || defined(_SPH_STR_SN)
     &                ,namestar
#endif

#ifndef _N_BODY
       real*8   u_ener,u_vol
#endif
!
!      global parameters:
!
!
      parameter(pi=3.141592654,two_third=.6666667)
#ifdef _TUBE
      parameter(xnorm=1.333333333)
#elif defined(_2D)
      parameter(xnorm=.4547284089)
#else
      parameter(xnorm=.3183098861)
#endif
      parameter(root_3=1.732050808)
      parameter (cc=0.744438)
      parameter (fine=0.56)
!
!            stars parameters:
#if !defined(_N_BODY) && ( defined(_SPH_STR) || defined(_SPH_STR_SN) )
      parameter (star_max=.044964,SN_mass=.036,solar=.0044964)
      parameter (star_min=.00144)
      parameter (u_sol=5.670e+01,u_ph=1.0e+02)
      parameter (r_sol=2.168e-08)
      parameter (SN_lif_time=3.0,SN_ene=1.573e+05)
#endif
      parameter(T1=3.37,T2=4.34,ucc=1.343295)
!
!
!      system array sizes:
!      -------------------
#if   defined(_1D)
      parameter(dim=1,vert=2)
      parameter(nn=4096,maxoctants=32760,nnode=4095)
#elif defined(_2D)
      parameter(dim=2,vert=4)
      parameter(nn=40960,maxoctants=163836,nnode=40959)
#else
      parameter(dim=3,vert=8)
# if   defined (_4096)
      parameter(nn=4096,maxoctants=32760,nnode=4095)
# elif defined (_8192)
      parameter(nn=8192,maxoctants=65528,nnode=8191)
# elif defined (_16384)
      parameter(nn=16384,maxoctants=131064,nnode=16383)
# elif defined (_16896)
      parameter(nn=16384,maxoctants=131064,nnode=16383)
# elif defined (_32768)
      parameter(nn=32768,maxoctants=262136,nnode=32767)
# elif defined (_40960)
      parameter(nn=40960,maxoctants=327672,nnode=40959)
# elif defined (_49152)
      parameter(nn=49152,maxoctants=393208,nnode=49151)
# elif defined (_65536)
      parameter(nn=65536,maxoctants=524280,nnode=65535)
# elif defined (_131072)
      parameter(nn=131072,maxoctants=1048568,nnode=131071)
# elif defined (_262144)
      parameter(nn=262144,maxoctants=2097144,nnode=262143)
# endif
#endif
!
!      usual integers and maximum number of neighbors:
!      -----------------------------------------------
      parameter(zero=0,one=1,two=2,three=3,four=4,eight=8)
#if   defined (_low_neighbor_resol)
      parameter(maxneighb=32)
#elif defined (_mid_neighbor_resol)
      parameter(maxneighb=96)
#elif defined (_high_neighbor_resol)
      parameter(maxneighb=128)
#elif defined (_very_high_neighbor_resol)
      parameter(maxneighb=256)
#endif

!      time resolution:
!      ----------------
#ifdef _time_resol_1
      parameter(ntbins=1)
#endif
#ifdef _time_resol_2
      parameter(ntbins=2)
#endif
#ifdef _time_resol_3
      parameter(ntbins=4)
#endif
#ifdef _time_resol_4
      parameter(ntbins=6)
#endif
#ifdef _time_resol_5
      parameter(ntbins=8)
#endif
#ifdef _time_resol_6
      parameter(ntbins=10)
#endif
#ifdef _time_resol_7
      parameter(ntbins=12)
#endif

#ifndef _N_BODY
CGS --> computer units:
      parameter(
     &      u_mass=  0.2248201E-35,
     &      u_leng=  0.3138732E-18,
     &      u_time=  0.3168568E-13,
     &      u_dens=  0.7270638E+20,
     &      u_den_n= 0.1207747E-03,
#ifndef _SPH_ADIABAT
     &      u_lum=   0.6962331E-32,
     &      u_cool=  0.2251604E+24,
#endif
     &      u_ener=  0.2206062E-45,
     &      u_vol=   0.3092165E-55
     &      )

Cloud abundances:
      parameter(XX=.75,YY=.25,ZZ=0)
#endif



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
#if !defined(_N_BODY) && ( defined(_SPH_STR) || defined(_SPH_STR_SN) )
     &                  ,namestar
#endif
!
!            PART II: 32-BITS SCALARS

#ifdef _DEBUG
      common /statement_counter/ sttmnt_count
#endif
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
#ifndef _N_BODY
     &                              ,dtm_safe
      common /fixed_neighbor_number/        n_fix
      common /neighb_guessing/        guessing
      common /neighb_toler/             sc,tol_nfix
      common /neighb_y/             y
      common /extras1/             alpha,beta,eta
#  if !defined(_N_BODY) && ( defined(_SPH_STR) || defined(_SPH_STR_SN) )
      common /star_number/             n_str,str_growth,n_fr
#  endif
#endif
#ifndef    _NON_SELFGRAVITATING
      common /extras2/theta,epsilon,theta2,eps2_max,eps2_min
#endif
      common /node/                   node,maxnode
      common /next/                   next
      common /total_number/           n
      common /total_mass/             tot_mass
      common /energies/              E_tot,E_min,E_max,Eo,W_tot
     &                              ,T_tot,U_tot,W_s
#ifdef _SPMHD
     &                              ,E_mag
#endif
#ifndef _N_BODY
     &                              ,W_p,luminosity
#endif
      common /Courant_Factor/             courant

#   ifdef _BACKGROUND_PRESSURE_
      common /bg_pressure/ Pbg
#   endif
!
!
!
!            PART III: SINGLE 3-D VECTORS
!
!
      common /rows/ row (two)
      common /momenta/         vmoment(dim)
     &                        ,angmom (dim)
#   ifdef _BACKGROUND_FIELD_
     &                        ,Bbg    (dim)
#   endif
      common /identity/       kronecker(dim,dim)
      common /levi/           levi_civita(dim,dim,dim)
      common /ifirst/         ifirst      (vert)
      common /count/          count      (vert)
!
!
!            PART IV: Rn-VECTORS
!
#ifndef _N_BODY
      common /hlist/            h_list(nn)
      common /status/       gas(nn)
#endif
      common /particle_mass/ mass(nn)
      common /particle_size/ s(nn)
#ifndef _N_BODY
      common /number_of_neighbours/       neighb            (nn)
      common /smoothing_lengths/        h            (nn)
     &                              ,h_dot            (nn)
      common /sph_densities/              rho            (nn)
      common /pressure/              P            (nn)
      common /temperatures/              Temp            (nn)
      common /speed_of_the_sound/        cs            (nn)
      common /thermal_energ/              
     &                               u            (nn)
     &                                ,u_o            (nn)
     &                              ,u_old            (nn)
     &                              ,u_dot            (nn)
     &                              ,u_dot_visc      (nn)
     &                              ,lambda            (nn)
     &                              ,vmonag            (nn)
#if !defined(_N_BODY) && ( defined(_SPH_STR) || defined(_SPH_STR_SN) )
      common /lifetime/              t_SN            (nn)
      common /star_masses/           star_mass      (nn)
#endif
      common /seeds/
     &                               iseed(nn)
     &                              ,jseed(nn)
     &                              ,kseed(nn)
     &                              ,lseed(nn)
     &                              ,mseed(nn)
     &                              ,nseed(nn)
      common /not_conv_list/          i_list(nn)
#endif
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
#ifndef _N_BODY
      common /sph_lists/       neighb_list(nn,maxneighb)
      common /viscotensor/       QQ(nn,maxneighb)
#      ifdef _SPH_PLUS_
      common /PLUSPLUS/       AddTerm(nn,maxneighb)
      common /PLUS/             d_grad_W(nn,maxneighb,dim)
#      endif
      common /grad_list/      grad_W(nn,maxneighb,dim)

#endif
!
!
!            PART VI: Rn3 VECTORS
!
      common /forces/
#if !(defined ( _N_BODY))
#      if ! defined (_TUBE) 
     &                         g    (nn,dim)
     &                        ,q_acc(nn,dim)
     &                        ,P_acc(nn,dim)
#      ifdef  _SPH_PLUS_
     &                        ,h_grad(nn,dim)
#      endif
#      else
     &                         q_acc(nn,dim)
     &                        ,P_acc(nn,dim)
#      endif
#else
     &                         g    (nn,dim)
#endif
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
#ifdef _COMPUTE_PHYSICAL_TIME
      common /potential_time/ time(ntbins)
#endif
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
#ifndef _TUBE
      common /gnodes/ gnode(maxoctants,dim)
#endif

      common /positions/ xp(nn,dim)
      common /velocities/ v(nn,dim)
#ifndef _N_BODY
      common /v0/ v_o(nn,dim), v_old(nn,dim)
      common /result_acc/ accel(nn,dim)
#endif

#ifdef _SPMHD
        common /B_field1/        B      (nn,dim)
        common /B_field2/        B_o      (nn,dim)
        common /B_field3/        XB      (nn,dim)
#endif

      common /boxy_coordinates/box      (nn,dim)

