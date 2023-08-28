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
     &            ,cc,accel,rho,v_o,v_old
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











      parameter(nn=40960,maxoctants=327672,nnode=40959)










!
!      usual integers and maximum number of neighbors:
!      -----------------------------------------------
      parameter(zero=0,one=1,two=2,three=3,four=4,eight=8)







      parameter(maxneighb=256)


!      time resolution:
!      ----------------







      parameter(ntbins=4)















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





      common /extras2/             theta,epsilon,theta2,epsilon2

      common /node/                   node,maxnode
      common /next/                   next
      common /total_number/           n
      common /total_mass/             tot_mass
      common /energies/              E_tot,E_min,E_max,Eo,W_tot
     &                              ,T_tot,U_tot,W_s

     &                              ,E_mag


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

      common /v0/ v_o(nn,dim), v_old(nn,dim)
      common /result_acc/ accel(nn,dim)



        common /B_field1/        B      (nn,dim)
        common /B_field2/        B_o      (nn,dim)
        common /B_field3/        XB      (nn,dim)


      common /boxy_coordinates/box      (nn,dim)

# 1 "Cool_interp.F"
      function cooling(n,T)
	implicit real (a-z)
	common/H_H2_CO/ Lamb(0:1023,0:63),Hy(0:1023,0:63)
	integer i_n,i_T
	if(n.eq.0)then
	  cooling=0
	  return
	endif
	if(T.eq.0)then
	  cooling=0
	  return
	endif
	lgn=log10(n)
	lgT=log10(T)
	!print*,lgn,lgT
	x_n=(lgn+2)*64/13
	x_T=lgT*1024/7
	!print*,x_n,x_T
	i_n=int(x_n)
	if(i_n.lt.0)then
	  cooling=0
	  return
	endif
	i_T=int(x_T)
	if(i_T.lt.0)then
	  cooling=0
	  return
	endif
	if(i_n.gt.63)i_n=63
	if(i_T.gt.1023)i_T=1023
	!print*,i_n,i_T
	L0=Lamb(i_T,i_n)
	i_n1=i_n
	if(i_n.lt.63)i_n1=i_n+1
	D_L_n=Lamb(i_T,i_n1)-L0
	i_T1=i_T
	if(i_T.lt.1023)i_T1=i_T+1
	D_L_T=Lamb(i_T1,i_n)-L0
	D_n=x_n-i_n
	D_T=x_T-i_T
	D_L=D_L_n*D_n+D_L_T*D_T
	L=L0+D_L
	L=10**L
	!print*,L
	!warning: positively returned value.
	! set efficience to be 25%
	!*** L=.25*L
	cooling=L
	return
      end
!
!
      subroutine Assembl_Table
	implicit real (a-z)
	common/H_H2_CO/ Lamb(0:1023,0:63),Hy(0:1023,0:63)
	integer i_n,i_T
	open(1,file='H+H2+CO.dat',status='old')
	print*,'Reading the cooling table.'
	print*,'Please, wait.'
	do i_n=0,63
	  do i_T=0,1023
	    read(1,*)x,x,Lamb(i_T,i_n),Hy(i_T,i_n)
	  end do
	end do
	close(1)
	print*,'Ready.'
	print*,'Writing the binary cooling table.'
	print*,'Please, wait.'
	open(1,file='H+H2+CO.bin'
     &	,status='unknown',form='unformatted')
	write(1)Lamb
	close(1)
	print*,'Ready.'
	return
      end
!
      function H_mol(n,T)
	implicit real (a-z)
	common/H_H2_CO/ Lamb(0:1023,0:63),Hy(0:1023,0:63)
	integer i_n,i_T
	if(n.eq.0)then
	  H_mol=0
	  return
	endif
	if(T.eq.0)then
	  H_mol=0
	  return
	endif
	lgn=log10(n)
	lgT=log10(T)
	x_n=(lgn+2)*64/13
	x_T=lgT*1024/7
	i_n=int(x_n)
	if(i_n.lt.0)then
	  H_mol=0
	  return
	endif
	i_T=int(x_T)
	if(i_T.lt.0)then
	  H_mol=0
	  return
	endif
	if(i_n.gt.63)i_n=63
	if(i_T.gt.1023)i_T=1023
	H0=Hy(i_T,i_n)
	i_n1=i_n
	if(i_n.lt.63)i_n1=i_n+1
	D_H_n=Hy(i_T,i_n1)-H0
	i_T1=i_T
	if(i_T.lt.1023)i_T1=i_T+1
	D_H_T=Hy(i_T1,i_n)-H0
	D_n=x_n-i_n
	D_T=x_T-i_T
	D_H=D_H_n*D_n+D_H_T*D_T
	H=H0+D_H
	H_mol=H
	return
      end
!
      subroutine Get_Table
	implicit real (a-z)
	common/H_H2_CO/ Lamb(0:1023,0:63),Hy(0:1023,0:63)
	open(1,file='H+H2+CO.bin'
     &	,status='old',form='unformatted')
	read(1)Lamb
	close(1)
	return
	end
!
      function heating(n_H,n_HI)
	implicit real (a-z)
	parameter(HR=3.8e-29)
	parameter(Hd=2.2e-28)
	parameter(zH=0.1)
	!print*,n_H,' cm^-3',n_HI,' cm^-3'
	Gamm_HR=HR*n_H
	Gamm_Hd=Hd*zH*n_H*n_HI
	Gamma= Gamm_HR + Gamm_Hd
	!print*,Gamma,' ergs cm^-3 sec^-1'
	!warning: positively returned value.
	!*** heating = Gamma
	heating = Gamma
	!*** heating = 4*Gamma
	!*** heating = 8*Gamma
	!*** heating = 32*Gamma
	!*** heating = 48*Gamma
	return
	end
# 1 "SPH_advance.F"


!macro definitions

      subroutine advance_energies
      INCLUDE 'SPH_common.h'
      parameter(el=.99375)

      print*,'% advance_energies'

      call Save_old_energies
      call SPH_first_law
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(gas(i).or.i_list(i))then
          u(i) = u_o(i) + .25*(dt_old(i)+dt(timebin))*u_dot(i)
          if(u(i).lt.0)then
            print*,
     &      'advance_energies: ERROR: negative thermal specific energy.'
            print*,'particle:',i
            print*,'Energy values:',u_o(i),u(i)
            print*,'time-bin:',timebin
            print*,'time-step:',dt(timebin)
            u(i)=u_o(i)
            if(u(i).lt.0)then
              print*,'advance_energies: ERROR:',
     &        ' cannot find a positive definite thermal energy'
              print*,'changing initial and final values.'
              u_o(i)=2*abs(u(i))
              u(i)=u_o(i)
            endif
          endif
        endif
      end do

      print*,'% advance_energies: done'

      end



# 75 "SPH_advance.F"


      subroutine Predict_velocities(factor,flag)
      integer flag
      INCLUDE 'SPH_common.h'
      if(flag.eq.1)then
          do i_=1,n_p(timebin)
          i=p_list(i_,timebin)
          if(gas(i).or.i_list(i))then
            do j=1,dim
              v(i,j) = .5*(factor*v(i,j)+(2-factor)*v_o(i,j))
            end do
          endif
        end do
      else
        call Hydrodyn_accel



        do i_=1,n_p(timebin)
          i=p_list(i_,timebin)
          if(gas(i))then
            do j=1,dim
              v(i,j)=v_o(i,j)+factor*.5*dt(timebin)*accel(i,j)
            end do
          endif
        end do
      endif
      end



      subroutine Move_particles
      INCLUDE 'SPH_common.h'
      do j=1,dim
        do i_=1,n_p(timebin)
          i=p_list(i_,timebin)
          dt_a =  .5*(dt(timebin)+dt_old(i))
          dt_b = .25*(dt(timebin)-dt_old(i))



          xp(i,j)=xp(i,j)+(v_o(i,j)+accel(i,j)*dt_b)*dt_a

        end do
      end do
      end



      subroutine advance_velocities
      INCLUDE 'SPH_common.h'

      print*,'advance_v'





      call save_old_velocities






      call Hydrodyn_accel









      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(i_list(i))then
          do j=1,dim
            v(i,j) = v_o(i,j) + accel(i,j) * dt(timebin)
          end do
        endif
      end do





      print*,'advance_v: done'

      end



      subroutine Correct_B
        INCLUDE 'SPH_common.h'
        do j=1,dim
          do i=1,n
            g(i,j)=v(i,j)
          end do
        end do
        call Predict_velocities(1.,1)
        call Advance_B
        do j=1,dim
          do i=1,n
            v(i,j)=g(i,j)
          end do
        end do
      end


      function B_rate(i,k)
      INCLUDE 'SPH_common.h'
      !---------------------------------------------------!
      !       Lagrangian form of the field evolution      !
      !---------------------------------------------------!
      !                                                   !
      !    $\dot{B}=(B\cdot\nabla){v}-(\nabla\cdot{V})B$  !
      !                                                   !
      ! where,                                            !
      !                                                   !
      !   B = b + B_0                                     !
      !                                                   !
      !---------------------------------------------------!

      divV_i=0
      Brate_i=0
      !for each neighbor of particle i do:
      do j_=1,neighb(i)
        j=neighb_list(i,j_)

        ! common scalar products:
        ! BdotNabla_ij = B^ dot grad_w^
        BdotNabla_ij=0
        ! VdotNabla_ij = v^ dot grad_w^
        VdotNabla_ij = 0
        do l=1,dim
          VdotNabla_ij = VdotNabla_ij + (v(j,l)-v(i,l)) * grad_w(i,j_,l)!Ok!

          BdotNabla_ij = BdotNabla_ij + (Bbg(l)+B(i,l)) * grad_w(i,j_,l)!Ok!



        end do !l

        divV_i = divV_i + VdotNabla_ij*mass(j) !Ok!
        Brate_i = Brate_i + BdotNabla_ij*mass(j)*(v(j,k)-v(i,k)) !Ok!
      end do !j

      ! returning result:

      B_rate =(Brate_i-(Bbg(k)+B(i,k))*divV_i)/rho(i)



      end


      subroutine Advance_B
      logical repeat
      INCLUDE 'SPH_common.h'
      parameter (tol=1e-5)

      print*,'Advance_B'

      do k=1,dim
        do i_=1,n_p(timebin)
          i=p_list(i_,timebin)
          if(gas(i))then

            B_o(i,k)=B(i,k)
            repeat=.true.
            do while(repeat)
              tau = .25 * (dt(timebin)+dt_old(i))
              B_old=B(i,k)
              dB = tau * B_rate(i,k)
              B(i,k) = .5 * ( B_o(i,k)+dB + B_old )
              control=abs(B(i,k)-B_old)-abs(B(i,k)+B_old)*tol
              repeat=(control.gt.0)
            end do

          endif
        end do
      end do

      print*,'Advance_B: done.'

      end


# 332 "SPH_advance.F"

# 1 "SPH_control.F"




      subroutine Settings
      INCLUDE 'SPH_common.h'
      iiseed=417253423




      n_fix = 90







      alpha = 0.5
      beta  = 1.0
      eta   = 0.1
      courant = 2.

      tol_nfix = 1

      tol_nfix = max(int( 4./sqrt(float(n_fix)) ),1)

      sc = n_fix**(.3333333)








      end

      subroutine TimeSet
      INCLUDE 'SPH_common.h'
      dtm=dt(0)










      do j=1,6


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


      print*,'Get_Table'
      call Get_Table

      dE = 0

!     perform the simulation now:
      do time_level=1,ninew

          global_time_level = global_time_level + 1

          rtime = global_time_level*dt(0)


          write(*,'(a6,i5,a2,f15.7,a5)')
     &            'Step',time_level,' (',rtime,' [t])'






          E_old = E_tot

          if(dt(0).eq.0.0)print*,'simula: warning: dt(0) is null'
          print*,'Integration_scheme (dE)'







          call Integration_scheme (0)
          if(dt(0).eq.0.0)then
            call finalize
            stop
          endif






          E_max=amax(E_max,E_tot)
          E_min=amin(E_tot,E_min)

!          integral parameters saving:
!          ---------------------------

!          extract prefix from output-file-names:

          l3        = index(percent,' ')

          if((global_time_level.eq.1).or.(Eo.eq.0))then
             open(2,file=relat,status='unknown',access='append')
             write(2,*)'file ',percent(1:l3)
             Eo=E_tot
             write(2,*)'Total energy = ',Eo
             write(2,*)'Total gravitational energy = ',W_tot
             E_max=Eo
             E_min=E_max
             close(2)
          endif
          ee = (E_min-E_max)/Eo

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

          open(2,file=evolution,status='unknown',access='append')
          write(2,'(x,14e16.8)')
     &             rtime                 ! 1
     &            ,W_tot+U_tot+T_tot     ! 2
     &            ,U_tot                 ! 3
     &            ,T_tot                 ! 4
     &            ,W_tot                 ! 5

     &            ,W_tot+W_s+W_p+2*T_tot ! 6
     &            ,luminosity            ! 7




     &            ,(vmoment(j),j=1,dim)  ! 8..10

     &            ,(angmom(j), j=1,dim)  ! 11..13


     &            ,E_mag                 ! 14

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

# 1 "SPH_cool.F"
! 
! 
! 



      subroutine cooling_rate
      real Lamb,n_H,n_H1
      INCLUDE 'SPH_common.h'
        common/H_H2_CO/ Lamb(0:1023,0:63),Hy(0:1023,0:63)
      freedom = 5
      do il=1,n_p(timebin)
        i = p_list(il,timebin)
        if (gas(i).and.i_list(i)) then
          Temp(i) = Temperature(u(i),rho(i),P(i),w_m,HI,freedom)

c         if(Temp(i).lt.(5.0))then
          do while(Temp(i).lt.(5.0))
            print*,'particle',i,' in timebin',timebin,':'
            print*,'warning: super-cooling'
            print*,'T=',Temp(i),' K'
            print*,'doubling specific thermal energy'
            u(i)=2*u(i)
            Temp(i) = Temperature(u(i),rho(i),P(i),w_m,HI,freedom)
          end do
c         endif

          ! get the mass fraction of H species in H cm^3:
          n_H1 = XX*rho(i)/u_den_n
          !** if(n_H1.lt.1)then
          !**   rho(i)=u_den_n/XX
          !** endif
          ! get the number of H species per cm^3:
          n_H = w_mol(1.,0.,HI)*n_H1
          ! get the number of H atoms per cm^3:
          n_H1 = HI * n_H1
          ! get the mean molecular number-density:
          !** n_tot=rho(i)*w_m/u_den_n
          ! get the CGS cooling (ergs cm^-3 sec^-1):
          LAMBDA(i) = cooling(n_H,Temp(i))
          ! add the CGS heating (ergs cm^-3 sec^-1):
          LAMBDA(i) = LAMBDA(i) - HEATING(n_H,n_H1)
          ! convert to computer units and change signal
          if ( rho(i) .lt. (LAMBDA(i)*1.0e-16) )LAMBDA(i) = 0
          LAMBDA(i) = LAMBDA(i) * u_cool
          LAMBDA(i) = - LAMBDA(i) / rho(i)
          ! Lambda is now the net heating per unit mass:
        endif
      end do

      if (you_can) call damp_cooling_rate

      end





















      subroutine damp_cooling_rate
      integer il,i
      real a,c
      INCLUDE 'SPH_common.h'

      print*,'* damp_cooling_rate'

      do il=1,n_p(timebin)
        i = p_list(il,timebin)
        if(gas(i).and.i_list(i))then
          c = -LAMBDA(i)
          if(c.gt.0)then
            a = .25 * u(i) / dt(timebin) + u_dot(i)
            if(c.gt.abs(a))then
              c = a/c
              LAMBDA(i) = - a / sqrt ( 1 + c*c )
            else
              a = c/a
              LAMBDA(i) = - c / sqrt ( 1 + a*a )
            endif
          endif
        endif
      end do

      print*,'* damp_cooling_rate: done'

      end


# 1 "SPH_grav.F"

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



      if (np(node,cell).eq.1)then
            ep2 = amax(2*s(label(node,cell)),5e-2)
            ep2 = ep2*ep2
      else
            ep2 = amax(s2,1e-3)
      endif
      ep2 = amin(ep2,epsilon2)

!      square softened distance:
      r2 = r2 + ep2
!      softened distance:
      r  = sqrt(r2)
!      cubic softened distance:
      r3 = r2 * r
!      calculate pole-vector:
      do j=1,dim
        ps(1,j)=ps(1,j)/r2
      end do
!      --------------------------------
!      perform monopole approximation
!      ------------------------------
      phi0 = - cellmass(node,cell) / r
!      --------------------------------
!      perform quadrupole approximation
!      --------------------------------
!      diagonal sumation:
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



      ep2 = s2
      !** commented on March 13, 2000
      !** if(np(node,cell).eq.1)then
        !** ep2 = amax(2*s(label(node,cell)),5e-2)
        !** ep2 = ep2*ep2
      !** else
        !** ep2 = amax(s2,1e-3)
      !** endif
      ep2 = amin(ep2,epsilon2)

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


# 416 "SPH_grav.F"

# 1 "SPH_integral.F"
! 
! 
! 
      subroutine Integral_quantities
      INCLUDE 'SPH_common.h'

      call Clear_integrals

      call Treepotential
      E_tot = W_tot




      call Kinetic_energy(1.)

      call Angular_momentum

      E_tot = E_tot + T_tot

      call Thermal_energy
      E_tot = E_tot + U_tot
      call Pressure_virial

      end


      subroutine Clear_integrals
      INCLUDE 'SPH_common.h'
      E_tot=0
      do j=1,dim
        angmom(j)=0
      end do
      end

      subroutine Angular_momentum
      INCLUDE 'SPH_common.h'
      call Clear_integrals
      do j=1,dim
        vmoment(j) = 0
!        for all system particles do:
        do i=1,n

          v_j = .5*(v(i,j)+v_o(i,j))



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



      subroutine Thermal_energy
      INCLUDE 'SPH_common.h'
!      clear total-thermal-energy accumulator:
      U_tot = 0
!      clear luminosity accumulator:
      luminosity = 0
        do i=1,n
        if(gas(i))then
          !calculate total thermal-energy
          U_tot = U_tot + mass(i) * u(i)
          if(u(i).eq.0)then
            print*,'Thermal_energy: error: u=0'
            stop
          endif
          luminosity = luminosity - mass(i) * lambda(i)
        endif
      end do

      print*,'Thermal_energy =',U_tot

      end

      subroutine Pressure_virial
        INCLUDE 'SPH_common.h'
        W_p = 0
        do i=1,n
          if(gas(i))then
            s_i = 0
            do j_=1,neighb(i)
              j = neighb_list(i,j_)
              GWX_ij = GWX(j_,i,j)*mass(j)
              s_i = s_i + GWX_ij*(P(i)+P(j)+QQ(i,j_))
            end do
          endif
          W_p = W_p + mass(i)*s_i
        end do
        W_p = -.5*W_p
      end


      subroutine Kinetic_energy(factor)
      !Why I did place the magnetic energy calculations
      !here don't ask me!
      INCLUDE 'SPH_common.h'
      ! clear kinetic-energy accumulator:
      T_tot = 0

      E_mag = 0

      do j=1,dim
        do i=1,n

          v_j = .5*((2-factor)*v(i,j)+factor*v_o(i,j))



          v2 = v_j * v_j
          T_tot = T_tot + .5*v2*mass(i)


          Bij_2=B(i,j)+Bbg(j)



          Bij_2=Bij_2*Bij_2
          E_mag=E_mag+Bij_2*mass(i)/rho(i)

        end do
      end do

      E_mag = E_mag/(8*pi)
      E_tot=E_tot+E_mag

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

      if(n_p(timebin).lt.1)return
      call initial_saving

      call move_particles



      call set_no_convergence
      call reset_time_list
      call predict_velocities(1.,0)
      call preliminaries


      call treegravity




      call advance_B
      call field_contrib


      call integrate_on_u_v

      call correct_B










      call save_old_timesteps

      print*,'predictor_corrector: done'

      end















      subroutine Integrate_on_u_v
      common/general_counter/ icounter

      logical repeat
      repeat= .true.

      print*,'Integrate_on_u_v'

      icounter=0

      do while (repeat)
        icounter=icounter+1
        write(*,'(a25,i3,a)')'predictor-corrector (pass',icounter,')'
        call Step_counter

        call predict_velocities(1.,1)
        call advance_energies



        call advance_velocities
        call check_convergence(repeat)
      end do






















      print*,'Integrate_on_u_v: done'

      end


      subroutine Set_i_list
      INCLUDE 'SPH_common.h'
      do i_=1,n_p(timebin)
        i_list(p_list(i_,timebin))=.true.
      end do
      end


      subroutine Check_convergence(repeat)
      logical repeat
      INCLUDE 'SPH_common.h'
      itol=17
      tol=1
      do i=1,itol
        tol=.5*tol
      end do
      repeat = .false.






      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        if(gas(i).and.i_list(i))then
          if(u(i).gt.0)then
            control=abs(u(i)-u_old(i))-tol*u(i)
            i_list(i)=(control.gt.0)
            control=(abs(u(i)-u_old(i))-tol*u_old(i))
            i_list(i)=i_list(i).or.(control.gt.0)
          endif
          do j=1,dim
            control=abs(v(i,j)-v_old(i,j))-tol*v(i,j)
            i_list(i)=i_list(i).or.(control.gt.0)
            control=abs(v(i,j)-v_old(i,j))-tol*v_old(i,j)
            i_list(i)=i_list(i).or.(control.gt.0)
          end do
          repeat = (repeat.or.i_list(i)).and.(iter.lt.1024)
        endif
      end do
      if(iter.eq.1024)then
        print*,'WARNING: Too many iterations on integrating (u,v).'
      endif

      end



      subroutine Preliminaries
      INCLUDE 'SPH_common.h'

      print*,'Preliminaries'










      call Maketree
      call Search_lengths

      call Densities
      call Get_kernel_grad

      print*,'end Preliminaries'

      end

      subroutine Step_counter
      INCLUDE 'SPH_common.h'
      iter = iter + 1
      end

      subroutine Set_no_convergence
      integer il,i
      INCLUDE 'SPH_common.h'
      iter = 0
      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        i_list(i)=.true.
      end do
      end

      subroutine Initial_saving
      integer il,i,j
      INCLUDE 'SPH_common.h'
      do j=1,dim
        do il=1,n_p(timebin)
          i=p_list(il,timebin)
          v_o(i,j)=v(i,j)
        end do
      end do

      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        if (gas(i)) then
          tau=.25*(dt(timebin)+dt_old(i))
          do j=1,dim
            v(i,j) = v_o(i,j) - .5*dt_old(i)*accel(i,j)
            B_o(i,j) = B(i,j) + tau*B_rate(i,j)
          end do
        endif
      end do


      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        if(gas(i))then
          tau=.25*(dt(timebin)+dt_old(i))
          u_o(i)=u(i)+tau*u_dot(i)! u_dot is the old value
        endif
      end do
      end

      subroutine Save_old_velocities
      integer il,i
      INCLUDE 'SPH_common.h'
      do il=1,n_p(timebin)
       i=p_list(il,timebin)
       do j=1,dim
         if(gas(i).or.i_list(i)) v_old(i,j) = v(i,j)
       end do
      end do
      end

      subroutine Save_old_energies
      integer il,i
      INCLUDE 'SPH_common.h'
      do il=1,n_p(timebin)
       i=p_list(il,timebin)
       if(gas(i).or.i_list(i)) u_old(i) = u(i)
      end do

      end



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

          i_list(i) = .true.
          gas(i)    = .true.
        end do
        call Preliminaries
        call Sound_speed
        call Pressures



        call Get_h_dot
        you_can  = (ichoice.eq.0)




        starting = .false.
        not_yet  = .false.
        first    = .false.
        if(dt(0).eq.0.0)return
      else
        starting = .true.
      endif
      if(ichoice.eq.1)then

        call SPH_First_law


        call Accelerations







        call Start_list
        do i=1,n
          dt_old(i) = 0
        end do
        ichoice=0
        you_can = .true.
        call predictor_corrector(0)
      endif
      call Start_list

      call Set_i_list


      print*,'initialize: done'

      end



      subroutine Integration_scheme (dE)
      INCLUDE 'SPH_common.h'

      print*,'Integration_scheme'


      call StartSeeds

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

      print*,'successful opening ',namecore

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
     &             epsilon      !--> softening length
      write(3)
     &             E_max
      write(3) 
     &             Eo
      write(3) 
     &             E_min

      write(3)
     &             n_fix
      write(3)
     &             alpha
      write(3)
     &             beta
      write(3)
     &             eta

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

      write(3)
     &            accel
      write(3)
     &            u
      write(3)
     &            u_dot
      write(3)
     &            vmonag
      write(3)
     &            gas








      write(3) 
     &            sc

      write(3)
     &            B


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
     &             epsilon      !--> softening length
      read(1)
     &             E_max
      read(1) 
     &             Eo
      read(1) 
     &             E_min

      read(1)
     &             n_fix
      read(1)
     &             alpha
      read(1)
     &             beta
      read(1)
     &             eta

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

      read(1)
     &            accel
      read(1)
     &            u
      read(1)
     &            u_dot
      read(1)
     &            vmonag
      read(1)
     &            gas








      read(1) 
     &             sc

      read(1)
     &             B


      close(1)
      epsilon2=epsilon*epsilon
      theta2=theta*theta
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
      print *,'   softening length (negative for default):'
      read(*,'(a132)')lntxt
      read(*,*)epsilon

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



      print*,'enter the background (homogeneous) field components ',
     &'(computer units):'
      read(*,'(a132)')lntxt
      read(*,*)(Bbg(i),i=1,dim)




      print*,'enter the background pressure (computer units):'
      read(*,'(a132)')lntxt
      read(*,*)Pbg






      print *
      print *,'parameters ok!'
      theta2=theta*theta

      if(epsilon.lt.0) call calc_epsilon
      epsilon2 = epsilon*epsilon

!     find a convenient time-setp:
      if(dt(0).lt.0)then

        dtm = .125*sqrt(epsilon2*epsilon/Tot_Mass)
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

     &   ,u(i)                               ! 8

     &   ,(B(i,j),j=1,dim)                   ! 9 - 11


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
      !gamma_ = -2*Tot_Mass*Tot_Mass*float(n)**(-.3333333)
      gamma_ = -0.5*Tot_Mass*Tot_Mass*float(n)**(-.3333333)
      cont = .true.
      epsilon = 0
      call maketree
      do while(cont)
        print*,'epsilon_old =',epsilon
        epsilon_old = epsilon
        w_old = W_tot
        call treepotential
        epsilon  = gamma_/(W_tot+w_old)
        cont = abs(epsilon-epsilon_old).gt.(epsilon/32)
      end do
      print*,'epsilon =',epsilon
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

      print *,'number of neighbors              : ',n_fix
      print *,'alpha                            : ',alpha
      print *,'beta                             : ',beta
      print *,'eta                              : ',eta
      print *
      print *



          print *,'      dissipative simulation!'




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

      write(4,*)'number of neighbors             : ',n_fix
      write(4,*)'alpha                           : ',alpha
      write(4,*)'beta                            : ',beta
      write(4,*)'eta                             : ',eta


      write(4,*)'background magnetic field       : ',(Bbg(i),i=1,dim)



      write(4,*)'background pressure             : ',Pbg




      write(4,*)'      non-adiabatic simulation'


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



      character*132 pre2




      character chrd*4
      INCLUDE 'SPH_common.h'
      l1  = index(nameout,']')+1
      l2  = index(nameout(l1:132),'.')-1
      l   = l1+l2
      ext = nameout(l:132)
      call numstring(chrd,mchr,ifile)



      pre2 = nameout(1:l-1) 

      open(2,file=pre2,status='unknown')
      print *,'       ',pre2
























      do i=1,n
        if(gas(i))then



          rho_CGS=rho(i)/u_dens
          call y_weight(rho_CGS*XX,Temp(i),HI)! --> HI
          





              write(2,'(17e16.8)')






     &         mass(i)            !  1
     &        ,(xp(i,j),j=1,dim)  !  2 - 4
     &        ,(v(i,j),j=1,dim)   !  5 - 7
     &        ,u(i)               !  8

     &        ,(B(i,j),j=1,dim)   !  9 - 11

     &        ,rho(i)             !  9 ! 12
     &        ,P(i)*rho(i)*rho(i) ! 10 ! 13
     &        ,-lambda(i)         ! 11 ! 14
     &        ,Temp(i)            ! 12 ! 15
     &        ,cs(i)              ! 13 ! 16

     &        ,HI                 ! 14 ! 17











        endif
      end do



      close(2)

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

          if(gas(i))then
             write(4,'(3e12.4)')
     &        (xp(i,j),j=1,dim)
*    &       ,(v(i,j),j=1,dim)
*    &       ,u(i)
*    &       ,P(i)*rho(i)*rho(i)
*    &       ,-lambda(i)
*    &       ,Temp(i)
*    &       ,rho(i)
          endif
      end do
      do i=1,n,n_show
          if(.not.gas(i))then
            write(4,'(3e12.4)')
     &       (xp(i,j),j=1,dim)
*    &      ,(v(i,j),j=1,dim)
*    &      ,u(i)
*    &      ,P(i)*rho(i)*rho(i)
*    &      ,-lambda(i)
*    &      ,1e+4
*    &      ,1e+4
          endif




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

# 1 "SPH_mhd.F"



c Function BB_grad(i,j,k) takes the k-component of the SPH-gradient from the
c Maxwell's stress-tensor.

      function BB_grad(i,jnb,j)
      INCLUDE 'SPH_common.h'

      
      B2=0
      do k=1,dim
        B2=B2+B(i,k)*B(i,k)
      end do

      
      Bscal=0
      do k=1,dim
        Bscal=Bscal+B(i,k)*grad_w(i,jnb,k)
      end do

      


      BB_grad=(B2*grad_w(i,jnb,j)-2*B(i,j)*Bscal)/
     &(8*pi*rho(i)*rho(i))

      end




      subroutine field_contrib
      










      
      INCLUDE 'SPH_common.h'
      do k=1,dim
        do i_=1,n_p(timebin)
          i=p_list(i_,timebin)
          if(gas(i))then
            XB(i,k)=0
            do jnb=1,neighb(i)
              j=neighb_list(i,jnb)
              XB(i,k)=XB(i,k) + mass(j)*(BB_grad(i,jnb,k)
     &                        + BB_grad(j,jnb,k))
            end do
            XB(i,k) = -XB(i,k)
          endif
        end do
      end do

      call background_field

      end


      subroutine background_field
      include 'SPH_common.h'
      do k=1,dim
        do i_=1,n_p(timebin)
          i=p_list(i_,timebin)
          if(gas(i))then
            contrib=0
            contrib1=0
            contrib2=0
            contrib3=0
            contrib4=0

            !scl2 -> Bbg^ dot B^_i/rho_i^2
            scl2=0
            do l=1,dim
              scl2=scl2+Bbg(l)*B(i,l)
            end do
            scl2=scl2/(rho(i)*rho(i))

            !for each neighbor do:
            do jnb=1,neighb(i)
              j=neighb_list(i,jnb)

              !common scalar products:
              !scl1 -> Bbg^ dot Grad^_i W_ij
              !scl3 -> Bbg^ dot B^_j
              scl1=0
              scl3=0
              do l=1,dim
                scl1=scl1+Bbg(l)*grad_w(i,jnb,l)
                scl3=scl3+Bbg(l)*B(j,l)
              end do! l

              !contrib1 -> Sum_j m_jBbg^dotGrad^_iW_ij B_j/rho_j^2(k-direction)
              contrib1=contrib1+
     &                  mass(j)*scl1*
     &                  B(j,k)/(rho(j)*rho(j))
              !contrib2 -> Sum_j m_j Bbg^dot Grad^_iW_ij
              contrib2=contrib2+mass(j)*scl1
              !contrib3 -> Sum_j m_j Grad^_iW_ij Bbg^dot B_j/rho_j^2
              contrib3=contrib3+
     &                  mass(j)*grad_w(i,jnb,k)*
     &                  scl3/(rho(j)*rho(j))
              !contrib4 -> Sum_j m_j Grad^_iW_ij (k-direction)
              contrib4=contrib4+
     &                  mass(j)*grad_w(i,jnb,k)
            end do! each neighbor

            contrib2=contrib2*B(i,k)/(rho(i)*rho(i))
            contrib4=contrib4*scl2
            contrib=(contrib1+contrib2-(contrib3+contrib4))/(4*pi)
            XB(i,k)=XB(i,k)+contrib
          endif
        end do
      end do
      end




# 1 "SPH_nns.F"

! 
!         SPH AND TREE-NEAREST-NEIGHBORS-SEARCHING SUBROUTINES
! 


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
      !!print*,'@ Set_permit'
      do node=1,maxnode
        permit(node)=.true.
      end do
      !!print*,'@ Set_permit: done'
      end


      subroutine NNS_gather(i)
      INCLUDE 'SPH_common.h'
      neighb(i) = 0
      call Set_permit
      do node=1,maxnode
        if(permit(node))then
          call Search_gather_neighbors(i,2*h(i))
        else
          call forbid(node)
        endif
      end do
      !!print*,i,neighb(i)
      end


      subroutine Search_gather_neighbors(i,r_i)
      




      integer cell
      INCLUDE 'SPH_common.h'
      r_cell = arest(node) + r_i
      do cell=1,vert
        if(np(node,cell).gt.0)then
          if (np(node,cell).eq.1) then
            if(spacing(cell,i,r_i,1).lt.0)neighb(i)=neighb(i)+1
          else
            if (spacing(cell,i,r_cell,0).gt.0) then
!             forbid tree descents by this path:
              next = down(node,cell)
              permit(next) = .false.
            endif
          endif
        endif
      end do
      return
      end



      subroutine neighboring
      INCLUDE 'SPH_common.h'
      !!print*,'! neighboring'
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        neighb(i) = 0
        call Set_permit
        do node=1,maxnode
          if(permit(node))then
            call effective_neighbors(i,2*h(i))
          else
            call Forbid(node)
          endif
        end do
      end do
      !!print*,'! neighboring: done'
      end



      subroutine effective_neighbors(i,r_i) 
      













      integer cell
      INCLUDE 'SPH_common.h'
      r_cell = arest(node) + r_i
      do cell=1,vert
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

      print*,'& Search_lengths'

      do i_=1,n_p(timebin)
        call smoothing_length( p_list(i_,timebin) )
      end do

      print*,'& Search_lengths: done'

      end


      subroutine smoothing_length(i)
      INCLUDE 'SPH_common.h'
!      initialize control variables:
      y           = 100
      guessing    = .3427301914159
      iterplus    = 0
      itercounter = 0
      converging  = .false.
!      start iterations counter:
      iter_neighb = 0
      !!print*,i,h(i),neighb(i)
      do while(.not.converging)
!        predict the smoothing-length:



        h(i)= h(i)*(1+y*((1+y)*n_fix/(n_fix+neighb(i)*y))**.3333333)/(1+y)

        !!print*,i,h(i),neighb(i)
        call NNS_gather(i)
        call Ask_for_convergence(i)
        iter_neighb = iter_neighb + 1
      end do
      !!print*,i,h(i),neighb(i)
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







# 438 "SPH_nns.F"




# 1 "SPH_sseeds.F"

      subroutine StartSeeds
!
      INCLUDE 'SPH_common.h'
!
      open(5,file='SPH_seeds.dat',status='old')
      read(5,*)iroot
      close(5)
      do i=1,n
            iseed(i) = 1791123513*ran(iroot)
            iseed(i)=iseed_max(iseed(i))
            seed=seed + iseed(i)/float(n)
            jseed(i) = 2008235129*ran(iroot)
            jseed(i)=iseed_max(jseed(i))
            seed=seed + jseed(i)/float(n)
            kseed(i) = 1235213459*ran(iroot)
            kseed(i)=iseed_max(kseed(i))
            seed=seed + kseed(i)/float(n)
            lseed(i) = 1611235341*ran(iroot)
            lseed(i)=iseed_max(lseed(i))
            seed=seed + lseed(i)/float(n)
            mseed(i) = 1352631411*ran(iroot)
            mseed(i)=iseed_max(mseed(i))
            seed=seed + mseed(i)/float(n)
            nseed(i) = 1251336413*ran(iroot)
            nseed(i)=iseed_max(nseed(i))
            seed=seed + nseed(i)/float(n)
      end do
      iroot=int(seed/6.)
      iroot=iseed_max(iroot)
      open(5,file='SPH_seeds.dat',status='unknown')
      !write(5,*)iroot
      close(5)
      return
      end

      function iseed_max(seed)
      integer seed
      seed=max(seed,2008235129-seed)
      seed=abs(float(seed))
      seed=seed+1-mod(seed,2)
      iseed_max=seed
      return
      end


# 1 "SPH_thermo.F"


! 
! 
! 


































      subroutine Pressures
      INCLUDE 'SPH_common.h'




      print*,'Pressures'











      do i_=1,n_p(timebin)
        i = p_list(i_,timebin)
        if(gas(i).and.i_list(i)) then

          !temperature --> T, P, 1/mu, HI, degrees of freedom
          Temp(i)=Temperature(u(i),rho(i),P(i),w_m,HI,freedom)





c         print*,u(i),rho(i)*P(i),w_m,HI,freedom
c         print*,'T=',Temp(i),i

        endif
      end do

      print*,'Pressures: done'

      end



      subroutine SPH_first_law
      INCLUDE 'SPH_common.h'



      call Cooling_rate !--> Temp, P, Lambda

      call Monaghan_tensor
      call Viscosity_heating
      call Adiabatic_heating

      call Total_rate

      end



      subroutine Total_rate
      INCLUDE 'SPH_common.h'
      do i_=1,n_p(timebin)
        i = p_list(i_,timebin)
        if(gas(i).and.i_list(i))u_dot(i) = lambda(i) + u_dot(i)
      end do
      end



      subroutine adiabatic_heating
      INCLUDE 'SPH_common.h'

      print*,'# adiabatic_heating'

      do i_=1,n_p(timebin)
        i = p_list(i_,timebin)
        if(gas(i).and.i_list(i))then



          u_left = 0
          u_rght = 0
          do j_=1,neighb(i)
            j = neighb_list(i,j_)
            s_left = 0
            s_rght = 0
            do k=1,dim
              s_left = s_left + grad_w(i,j_,k) * v(i,k)
              s_rght = s_rght + grad_w(i,j_,k) * v(j,k)
            end do
            GWV_ij = (s_left-s_rght) * mass(j)
            u_left = u_left + GWV_ij
            u_rght = u_rght + GWV_ij * P(j)




          end do
          u_left = u_left * P(i)
          u_dot(i) = .5*(u_left + u_rght)+u_dot_visc(i)



        endif
      end do

      print*,'# adiabatic_heating: done'

      end


      subroutine viscosity_heating
      INCLUDE 'SPH_common.h'

      print*,'# viscosity_heating'

      do i_=1,n_p(timebin)
        i = p_list(i_,timebin)
        if(gas(i).and.i_list(i))then
          u_dot_visc(i) = 0
          do j_=1,neighb(i)
            j = neighb_list(i,j_)
            s_left = 0
            s_rght = 0
            do k=1,dim
              s_left = s_left + grad_w(i,j_,k) * v(i,k)
              s_rght = s_rght + grad_w(i,j_,k) * v(j,k)
            end do
            GWV_ij = ((s_left-s_rght)*QQ(i,j_)) * mass(j) 
            u_dot_visc(i) = u_dot_visc(i) + GWV_ij
          end do
          u_dot_visc(i) = .5 * u_dot_visc(i)
        endif
      end do

      print*,'# viscosity_heating: done'

      end




      function GWX(j_,i,j)
      INCLUDE 'SPH_common.h'
      s_left = 0
      s_rght = 0
      do k=1,dim
        s_left = s_left + grad_w(i,j_,k) * xp(i,k)
        s_rght = s_rght + grad_w(i,j_,k) * xp(j,k)
      end do
      GWX = s_left-s_rght
      end

      subroutine Sound_speed




!     this procedure calculates the adiabatic speed of sound:

            ! gamma - 1 = ( P/rho )/u
            ! c^2 = gamma ( gamma - 1 ) u
            !     = (P/rho) [1 + P/(rho u)]

      INCLUDE 'SPH_common.h'

      print*,'# sound'

      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(gas(i).and.i_list(i))then
          p_by_rho=P(i)*rho(i)

          


          ! cs(i) = u(i)+p_by_rho
          


          cs(i) = 8.154319e-3     *(3*(0.75            +.5*0.25            )*Temp(i)+5.1955e+4       *0.75            )*0.5








          if (cs(i).lt.0) then
            print*,flsh,'error: invalid speed of sound!',norm
            print*,'c^2 =',cs(i)
            stop
          endif

          







          





          




          v2=0
          do k=1,dim
            v2=v2+v(i,k)*v(i,k)
          end do
          B2=0
          do k=1,dim
            Bk=0
            do j=1,dim
            do l=1,dim
              
              Bk=Bk+levi_civita(k,j,l)*B(i,j)*v(i,l)
            end do
            end do
            B2=B2+Bk*Bk
          end do
          B2=B2/(v2+0.0005)
          cs(i) = cs(i) + B2/(8*pi*rho(i)+0.0005)

          cs(i)=sqrt(cs(i))
        endif
      end do

      print*,'# sound: done'

      end



      subroutine Monaghan_tensor
      INCLUDE 'SPH_common.h'
      call Sound_speed

      print*,'# Monghan Tensor'

      do i_=1,n_p(timebin)
        i = p_list(i_,timebin)
	!print*,'i=',i,'i_=',i_
        if(gas(i).and.i_list(i))then
          vmonag(i) = 0
          do j_=1,neighb(i)
            j=neighb_list(i,j_)
	    !print*,'j=',j,'j_=',j_
            Vmu_ij=0
            do k=1,dim
              x_ij=xp(i,k)-xp(j,k)
              v_ij=v(i,k)-v(j,k)
              Vmu_ij=Vmu_ij+x_ij*v_ij
              !print*,x_ij,v_ij,vmu_ij
            end do
            !print*,vmu_ij
            if(vmu_ij.lt.0)then
              h_ij = .5 * ( h(i) + h(j) )
	      !print*,'h_ij=',h_ij
              h_eta_2 = eta * h_ij
              h_eta_2 = h_eta_2 * h_eta_2
              r_2 = 0
              do k=1,dim
                x_ij = xp(i,k) - xp(j,k)
                r_2  = r_2  + x_ij*x_ij
              end do
              over_r_2 = 1 / ( h_eta_2 + r_2 )
              Vmu_ij = - Vmu_ij * h_ij * over_r_2
              c_ij=.5*(cs(i)+cs(j))
              rho_ij=.5*(rho(i)+rho(j))
              Vmu_ij = (alpha*c_ij + beta*Vmu_ij)*Vmu_ij
              vmonag(i)=amax(vmonag(i),Vmu_ij)
              QQ(i,j_) = Vmu_ij/rho_ij
            else
              QQ(i,j_) = 0
            endif
          end do
        endif
      end do

      print*,'# Monaghan tensor: done.'

      end

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

      B2=0

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

      if(first)then























        do i=1,n
          h(i)=0.3*s(i)
        end do

          first = .false.
      endif


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


!
!

























# 1 "SPH_usr.F"

! 
! 
! 









      subroutine Densities
      




      INCLUDE 'SPH_common.h'

      print*,'|-> densities'

      call Neighboring
      do i_=1,n_p(timebin)
        if(gas(p_list(i_,timebin)))call Rho_(p_list(i_,timebin))
      end do

      print*,'|-> densities: done'

      end



      subroutine Rho_(i)
      INCLUDE 'SPH_common.h'
      rho(i)=0
      do j_=1,neighb(i)
        j=neighb_list(i,j_)
        r_ij = 0
        do k=1,dim
          x_ij = xp(i,k)-xp(j,k)
          r_ij = r_ij + x_ij*x_ij
        end do
        r_ij = sqrt(r_ij)
        
        w_ij = .5*(sphKernel(r_ij,h(i))+sphKernel(r_ij,h(j)))
        rho(i) = w_ij * mass(j) + rho(i)
      end do

      



      if(rho(i).lt.(1.2e-4))rho(i)=1.2e-4










      end



      subroutine Get_kernel_grad
      INCLUDE 'SPH_common.h'
      common /aux/ aux(nn,dim)
      equivalence (aux,x_ij)
      dimension x_ij(nn,dim)
      

      print*,'kernel gradients'

      do i_=1,n_p(timebin)
        i = p_list(i_,timebin)
        if(gas(i))then
          do j_=1,neighb(i)
            j=neighb_list(i,j_)
            r_ij = 0
            do k=1,dim
              x_ij(1,k) = xp(i,k) - xp(j,k)
              r_ij = r_ij + x_ij(1,k)*x_ij(1,k)
            end do
            r_ij=sqrt(r_ij)
            if(r_ij.gt.0)then
            





              grad_ww=.5*(D_kernel(r_ij,h(i))+D_kernel(r_ij,h(j)))
            



              do k=1,dim
                grad_w(i,j_,k) = x_ij(1,k) * grad_ww / r_ij
              end do
            else
              do k=1,dim
                grad_w(i,j_,k) = 0
              end do
            endif
          end do
        endif
      end do

      print*,'kernel gradients: done!'





      end






















      subroutine Get_h_dot
      



      INCLUDE 'SPH_common.h'
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(h_list(i).and.gas(i))then
          
          h_dot(i) = 0
          do j_=1,neighb(i)
            j=neighb_list(i,j_)
            if (gas(j)) then
              do k=1,dim
                h_dot(i)=h_dot(i)+(v(j,k)-v(i,k))*mass(j)*grad_w(i,j_,k)
              end do
            endif
          end do
          
          tau=.25*(dt(timebin)+dt_old(i))
          
          h_dot(i) = h_dot(i) * tau / (rho(i) * 3)
        endif
      end do
      end




# 247 "SPH_usr.F"


# 1 "SPH_usr_accel.F"

      subroutine Pressure_accel
      INCLUDE 'SPH_common.h'
      do k=1,dim
        do i_=1,n_p(timebin)
          i = p_list(i_,timebin)
          if(gas(i).and.i_list(i))then
            P_acc(i,k) = 0
            P_k = 0
            q_acc(i,k) = 0
            do j_=1,neighb(i)
              j=neighb_list(i,j_)
              q_acc(i,k) = mass(j)*(grad_w(i,j_,k)*QQ(i,j_)) + q_acc(i,k)
              gw_ijk     = mass(j)*grad_w(i,j_,k)


              P_k        = P_k + gw_ijk*(P(j)-Pbg/(rho(j)*rho(j)))






              P_acc(i,k) = P_acc(i,k) + gw_ijk
            end do
            q_acc(i,k) = -q_acc(i,k)


            P_acc(i,k) = -(P_k+P_acc(i,k)*(P(i)-Pbg/(rho(i)*rho(i))))






          endif
        end do
      end do





















      end


      subroutine Accelerations 


      call treegravity




      call Hydrodyn_accel
      end



      subroutine Hydrodyn_accel
      call Pressure_accel

      call Field_contrib

      call Result_accel
      end

 

      subroutine Result_accel
      INCLUDE 'SPH_common.h'
      do j=1,dim
        do i_=1,n_p(timebin)
          i = p_list(i_,timebin)
          if (i_list(i).and.gas(i)) then
                   accel(i,j) = q_acc(i,j) + P_acc(i,j)


     &               + XB(i,j)

            accel(i,j) = accel(i,j) + g(i,j)

          endif
        end do
      end do



      end

