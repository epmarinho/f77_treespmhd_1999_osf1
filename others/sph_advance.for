# 1 "SPH_advance.F"


      subroutine advance_energies
      INCLUDE 'SPH_common.h'
      parameter(el=.99375)
      call Save_old_energies
      call SPH_first_law
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(gas(i).or.i_list(i))then
          u(i) = u_o(i) + .25*(dt_old(i)+dt(timebin))*u_dot(i)
          if(u(i).lt.(0.0625))then
            u_o(i)=.0625
            u(i)=u_o(i)
          endif
        endif
      end do
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
      call save_old_velocities
      call Pressure_accel
      call Result_accel
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(i_list(i))then
          do j=1,dim
            v(i,j) = v_o(i,j) + accel(i,j) * dt(timebin)
          end do
        endif
      end do
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
      !* print*,k,i,Brate_i
      end


      subroutine Advance_B
      logical repeat
      INCLUDE 'SPH_common.h'
      parameter (tol=1e-5)
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
      end


# 326 "SPH_advance.F"
