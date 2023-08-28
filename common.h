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
     &             ,time_level,global_time_level,timebin,p_list,label,vert
     &             ,deepest,nn,maxoctants,nnode,maxneighb,ichoice,is,ns
     &             ,lst,ifile,ni
     &             ,kwsc,n_timebins,ntbins2,n,next,maxnode,n_p,np,ifirst,ini,nit
     &             ,iiseed,kronecker,levi_civita,sttmnt_count

     &             ,iter,itercounter,iterplus,iter_neighb
     &             ,n_fix,neighb_list,neighb
     &             ,iseed,jseed,kseed,lseed,mseed,nseed,tol_nfix




      real     mass,cellmass,epsilon,g,arest,pi,epsilon2,dt_elapsed,dt,dtm
     &            ,dt_old,p_inertia,x_cell
     &            ,angmom,s,v,xp,theta,theta2,tot_mass
     &            ,E_tot,W_tot,T_tot,U_tot,E_min,E_max,Eo,courant,vmoment,gnode

     &            ,lambda,luminosity,two_third,xnorm,root_3,sc,W_p,W_s
     &            ,cc,accel,rho,v_o
     &            ,u,u_dot,u_o,h,h_dot
     &            ,T1,T2,ucc,u_mass,u_leng
     &            ,fine,u_time,u_dens,u_den_n,XX,YY,ZZ,dtm_safe
     &            ,alpha,beta,eta,guessing
     &            ,P,QQ,P_acc,q_acc,grad_W,u_dot_visc,Temp,cs,u_old,vmonag

     &            ,u_lum,u_cool,u_ph,E_n


     &            ,r_sol,SN_lif_time,SN_ene,t_SN,star_mass,n_str,y,str_growth
     &            ,star_min,star_max,SN_mass,solar,u_sol


      logical  nextnode,permit,convergence
     &            ,you_can,starting
     &            ,converging,ok,nonzero,nonover,nz,nv,ovc
     &            ,not_yet,first

     &            ,h_list,i_list,gas

      character esc,cr,cutl*4,rev*4,flsh*4,hltd*4,norm*4,home*4,cls*4
     &              ,answer
      character*132 evolution,name,nameout,namecore,relat,percent





       real*8   u_ener,u_vol

!
!      global parameters:
!
!
      parameter(pi=3.141592654,two_third=.6666667)



      parameter(xnorm=.4547284089)



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




      parameter(dim=2,vert=4)
      parameter(nn=40960,maxoctants=163836,nnode=40959)
























!
!      usual integers and maximum number of neighbors:
!      -----------------------------------------------
      parameter(zero=0,one=1,two=2,three=3,four=4,eight=8)



      parameter(maxneighb=96)






!      time resolution:
!      ----------------



















      parameter(ntbins=12)



CGS --> computer units:
      parameter(
     &      u_mass=  0.2248201E-35,
     &      u_leng=  0.3138732E-18,
     &      u_time=  0.3168568E-13,
     &      u_dens=  0.7270638E+20,
     &      u_den_n= 0.1207747E-03,

     &      u_lum=   0.6962331E-32,
     &      u_cool=  0.2251604E+24,

     &      u_ener=  0.2206062E-45,
     &      u_vol=   0.3092165E-55
     &      )

Cloud abundances:
      parameter(XX=.75,YY=.25,ZZ=0)




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


      common /statement_counter/ sttmnt_count

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

     &                              ,dtm_safe
      common /fixed_neighbor_number/        n_fix
      common /neighb_guessing/        guessing
      common /neighb_toler/             sc,tol_nfix
      common /neighb_y/             y
      common /extras1/             alpha,beta,eta







      common /node/                   node,maxnode
      common /next/                   next
      common /n/                   n
      common /total_mass/             tot_mass
      common /energies/              E_tot,E_min,E_max,Eo,W_tot
     &                              ,T_tot,U_tot,W_s




     &                              ,W_p,luminosity

      common /Courant_Factor/             courant


      common /bg_pressure/ Pbg

!
!
!
!            PART III: SINGLE 3-D VECTORS
!
!
      common /rows/ row (two)
      common /momenta/         vmoment(dim)
     &                        ,angmom (dim)

     &                        ,Bbg    (dim)

      common /identity/       kronecker(dim,dim)
      common /levi/           levi_civita(dim,dim,dim)
      common /ifirst/         ifirst      (vert)
      common /count/          count      (vert)
!
!
!            PART IV: Rn-VECTORS
!

      common /hlist/            h_list(nn)
      common /status/       gas(nn)

      common /particle_mass/ mass(nn)
      common /particle_size/ s(nn)

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




      common /seeds/
     &                               iseed(nn)
     &                              ,jseed(nn)
     &                              ,kseed(nn)
     &                              ,lseed(nn)
     &                              ,mseed(nn)
     &                              ,nseed(nn)
      common /not_conv_list/          i_list(nn)

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
!        /  row  /|           j
!      -------------------------------------
!            i  |  1  |  2  |   3  |   4   
!      -------------------------------------
!            1  |  t  |label| cell | iflg  
!            2  |  "  | "   |  "   |   "   
!            3  |  "  | "   |  "   |   "   
!             ... | ... |...  | ...  |  ...  
!              n  |  "  | "   |  "   |   "   
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

      common /sph_lists/       neighb_list(nn,maxneighb)
      common /viscotensor/       QQ(nn,maxneighb)




      common /grad_list/      grad_W(nn,maxneighb,dim)


!
!
!            PART VI: Rn3 VECTORS
!
      common /forces/


     &                         g    (nn,dim)
     &                        ,q_acc(nn,dim)
     &                        ,P_acc(nn,dim)










!
!
!            PART VII: Rn nt VECTORS
!
!
      common /time_bin/             n_p(ntbins)
      common /timelist/             p_list(nn,ntbins)
      common /individual_time/        dt_old(nn)
      common /time_steps/              dt(0:ntbins)
      common /elapsed_time/              dt_elapsed(ntbins)
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


      common /positions/           xp      (nn,dim)
      common /velocities/          v        (nn,dim)

      common /v0/                   v_o      (nn,dim)
      common /result_acc/          accel      (nn,dim)








      common /boxy_coordinates/box      (nn,dim)

