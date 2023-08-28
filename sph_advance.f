#ifndef _N_BODY

!macro definitions
#      ifndef _SPH_ISOTHERM
      subroutine advance_energies
      INCLUDE 'SPH_common.h'
      parameter(el=.99375)
#ifdef _VERBOSE
      print*,'% advance_energies'
#endif
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
#ifdef _VERBOSE
      print*,'% advance_energies: done'
#endif
      end
#        endif


#      ifdef _FAST_SEARCHING
      subroutine Advance_lengths /* Velocity sensitive */
      INCLUDE 'SPH_common.h'
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        h_list(i)=.true.
        if ((h(i).gt.0).and.(abs(h_dot(i)).lt.0.5)) then
          h(i)=(1+h_dot(i))*h(i)
        else
          h_list(i)=.false.
        endif
      end do
      call Get_h_dot
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if (h_list(i).and.((h(i).gt.0).and.(abs(h_dot(i)).lt.0.5))) then
          h(i)=h(i)/(1-h_dot(i))
        else
          h(i) = h(i)*(1-h_dot(i))/(1+h_dot(i))
          if(h(i).le.0) h(i)=2*s(i)
#ifdef _VERBOSE
          print*,'      Redoing smoothing-length',h(i),' for particle',i
          print*,'      with',neighb(i),' neighbours'
#endif
          call smoothing_length(i)
#ifdef _VERBOSE
          print*,'      New smoothing-length',h(i),' for particle',i
          print*,'      with',neighb(i),' neighbours'
          print*,'      Done.'
#endif
        endif
      end do
      end
#      endif

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
#       ifdef _ROTATING_FRAME
        call Coriolis
#       endif _ROTATING_FRAME
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
#endif /* _N_BODY */


      subroutine Move_particles
      INCLUDE 'SPH_common.h'
      do j=1,dim
        do i_=1,n_p(timebin)
          i=p_list(i_,timebin)
          dt_a =  .5*(dt(timebin)+dt_old(i))
          dt_b = .25*(dt(timebin)-dt_old(i))
#ifdef _N_BODY
          xp(i,j)=xp(i,j)+(v(i,j)+g(i,j)*dt_b)*dt_a
#else
          xp(i,j)=xp(i,j)+(v_o(i,j)+accel(i,j)*dt_b)*dt_a
#endif
        end do
      end do
      end



      subroutine advance_velocities
      INCLUDE 'SPH_common.h'
#ifdef _VERBOSE
      print*,'advance_v'
#endif
#ifdef _N_BODY
      call maketree
      call treegravity
#else
      call save_old_velocities
#  ifdef _SPH_ISOTHERM
      call Pressures
#  endif
#  ifdef _VRM_INCOMPRESSIBLE_SPH
      call VRM
#  endif
      call Hydrodyn_accel
#endif
#ifdef _N_BODY
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        do j=1,dim
          v(i,j) = v(i,j) + g(i,j) * dt(timebin)
        end do
      end do
#else
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(i_list(i))then
          do j=1,dim
            v(i,j) = v_o(i,j) + accel(i,j) * dt(timebin)
          end do
        endif
      end do
#  ifdef _VRM_INCOMPRESSIBLE_SPH
      call VRM
#  endif
#endif
#ifdef _VERBOSE
      print*,'advance_v: done'
#endif
      end


#ifdef _SPMHD
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
#ifdef _BACKGROUND_FIELD_
          BdotNabla_ij = BdotNabla_ij + (Bbg(l)+B(i,l)) * grad_w(i,j_,l)!Ok!
#else
          BdotNabla_ij = BdotNabla_ij + B(i,l)*grad_w(i,j_,l) !Ok!
#endif
        end do !l

        divV_i = divV_i + VdotNabla_ij*mass(j) !Ok!
        Brate_i = Brate_i + BdotNabla_ij*mass(j)*(v(j,k)-v(i,k)) !Ok!
      end do !j

      ! returning result:
#ifdef _BACKGROUND_FIELD_
      B_rate =(Brate_i-(Bbg(k)+B(i,k))*divV_i)/rho(i)
#else
      B_rate =(Brate_i-B(i,k)*divV_i)/rho(i)
#endif
      end


      subroutine Advance_B
      logical repeat
      INCLUDE 'SPH_common.h'
      parameter (tol=1e-5)
#ifdef _VERBOSE
      print*,'Advance_B'
#endif
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
#ifdef _VERBOSE
      print*,'Advance_B: done.'
#endif
      end
#endif

#ifdef _VRM_INCOMPRESSIBLE_SPH
      subroutine VRM
c	This subroutine performs the Marinho's Velocity Relaxation Method
c	of simulating SPH incompressible fluid dynamics (Marinho 2000, in
c	preparation.
c
      INCLUDE 'SPH_common.h'
      parameter (tol=1e-5)
      !** !* print*,'VRM: start.'
      do iteration=1,5
      !foreach timebin-particle do
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(i_list(i))then
          !perform velocity relaxation in the neighborings:
          do j_=1,neighb(i)-1
            j=neighb_list(i,j_)
            do k_=j_+1,neighb(i)
              k=neighb_list(i,k_)
              eta_jk=0
              r2_jmk=0
              do l=1,dim
                r2_jk=r2_jk+(xp(j,l)-xp(k,l))*(xp(j,l)-xp(k,l))
                eta_jk=eta_jk+
     &          (v_o(j,l)-v_o(k,l))*
     &          (xp(j,l)-xp(k,l))
              end do !l
              eta_jk=eta_jk/(r2_jk*(mass(j)+(mass(k))))
              !to be optimized in future:
              do l=1,dim
                v_o(j,l)=v_o(j,l)-mass(k)*(xp(j,l)-xp(k,l))*eta_jk
                v_o(k,l)=v_o(k,l)+mass(j)*(xp(j,l)-xp(k,l))*eta_jk
              end do !l
            end do !k_
          end do !j_
        endif
      end do !i_
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(i_list(i))then
          !perform velocity relaxation in the neighborings:
          do j_=1,neighb(i)-1
            j=neighb_list(i,j_)
            do k_=j_+1,neighb(i)
              k=neighb_list(i,k_)
              eta_jk=0
              r2_jmk=0
              do l=1,dim
                r2_jk=r2_jk+(xp(j,l)-xp(k,l))*(xp(j,l)-xp(k,l))
                eta_jk=eta_jk+
     &          (v(j,l)-v(k,l))*
     &          (xp(j,l)-xp(k,l))
              end do !l
              eta_jk=eta_jk/(r2_jk*(mass(j)+(mass(k))))
              !to be optimized in future:
              do l=1,dim
                v(j,l)=v(j,l)-mass(k)*(xp(j,l)-xp(k,l))*eta_jk
                v(k,l)=v(k,l)+mass(j)*(xp(j,l)-xp(k,l))*eta_jk
              end do !l
            end do !k_
          end do !j_
        endif
      end do !i_
      end do ! iteration
      !** !* print*,'VRM: done.'
      end
#endif /* _VRM_INCOMPRESSIBLE_SPH */
