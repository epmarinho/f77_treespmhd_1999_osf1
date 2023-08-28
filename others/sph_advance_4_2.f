#ifndef _N_BODY
/* Local macro definitions */

#      ifndef _SPH_ISOTHERM
      subroutine Advance_energies
      INCLUDE 'SPH_common.h'
      parameter(el=.99375)
      call Save_old_energies
      call SPH_first_law
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(gas(i).or.i_list(i))then
          u(i) = u_o(i) + .25*(dt_old(i)+dt(timebin))*u_dot(i)
#                if defined (_SPH_ADIABAT)
          if(u(i).lt.0)then
            print*,
     &            'Advance_energies: ERROR: negative thermal specific energy.'
            print*,'particle:',i
            print*,'Energy values:',u_o(i),u(i)
            print*,'time-bin:',timebin
            print*,'time-step:',dt(timebin)
            u(i)=u_o(i)
            if(u(i).lt.0)then
              print*,'Advance_energies: ERROR:',
     &                   ' cannot find a positive definite thermal energy'
              /* stop */
              print*,'changing initial and final values.'
              u_o(i)=2*abs(u(i))
              u(i)=u_o(i)
                  endif
          endif
#                else
          /* check energy integrity: */
          if(u(i).lt.(0.0625))then
            u_o(i)=.0625
            u(i)=u_o(i)
          endif
#                endif
        endif
      end do
      end
#        endif


      subroutine Advance_lengths /* Velocity sensitive */
#      ifdef _FAST_SEARCHING
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
          print*,'      Redoing smoothing-length',h(i),' for particle',i
          print*,'      with',neighb(i),' neighbours'
          call particle_length(i)
          print*,'      New smoothing-length',h(i),' for particle',i
          print*,'      with',neighb(i),' neighbours'
          print*,'      Done.'
        endif
      end do
#      endif
      end

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



      subroutine Advance_v
      INCLUDE 'SPH_common.h'
#ifdef _N_BODY
      call maketree
      call treegravity
#else
#	ifdef _SPH_ISOTHERM
      call Pressures
#	endif
#	ifdef _VRM_INCOMPRESSIBLE_SPH
      call VRM
#	endif
      call Pressure_accel
      call Result_accel
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
#	ifdef _VRM_INCOMPRESSIBLE_SPH
      call VRM
#	endif
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
      !---------------------------------------------!
      !    Lagrangian form of the field evolution   !
      !---------------------------------------------------!
      !                                                   !
      !    $\dot{B}=(B\cdot\nabla){v}-(\nabla\cdot{V})B$  !
      !                                                   !
      ! where,                                            !
      !                                                   !
      !   B = b + B_0                                     !
      !                                                   !
      !---------------------------------------------------!
      div=0
      !** div2=0
      B_rate_k=0
      do j_=1,neighb(i)
        j=neighb_list(i,j_)
        ! Bscal = B^ dot grad_w^
        Bscal=0
        vscal = 0
        !** scal2= 0
        do l=1,dim
          vscal = vscal + (v(j,l)-v(i,l)) * grad_w(i,j_,l) !Ok!
          !** scal2= scal2 + (B(j,l)-B(i,l)) * grad_w(i,j_,l) !Ok!
          Bscal = Bscal + B(i,l)*grad_w(i,j_,l) !Ok!
        end do
        div      = div      + vscal*mass(j) !Ok!
        !** div2   =div2     + scal2*mass(j) !Ok!
        B_rate_k = B_rate_k+(Bscal*mass(j))*(v(j,k)-v(i,k)) !Ok!
      end do
      ! output result:
C     B_rate_k= (B_rate_k - B(i,k)*div + v(i,j)*div2 )/
C    &           (rho(i)+0.0000005)
      B_rate_k=(B_rate_k-B(i,k)*div)/(rho(i)+0.0000005)
      B_rate = B_rate_k
      end


      subroutine Advance_B
      logical repeat
      INCLUDE 'SPH_common.h'
      parameter (tol=1e-5)
      !print*,'Advance_B'
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
      !print*,'Done.'
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
      print*,'VRM: start.'
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
      print*,'VRM: done.'
      end
#endif /* _VRM_INCOMPRESSIBLE_SPH */
