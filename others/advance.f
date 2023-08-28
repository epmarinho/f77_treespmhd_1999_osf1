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
        call Coriolis
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



      subroutine advance_v
      INCLUDE 'SPH_common.h'
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
