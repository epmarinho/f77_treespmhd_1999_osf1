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
              P_k        = P_k + gw_ijk*(P(j)-Pbg/rho(j))
              P_acc(i,k) = P_acc(i,k) + gw_ijk
            end do
            q_acc(i,k) = -q_acc(i,k)
            P_acc(i,k) = -(P_k+P_acc(i,k)*(P(i)-Pbg/rho(i)))
          endif
        end do
      end do
      end

      subroutine Accelerations 
      call extern_gravity_2d
      call Hydrodyn_accel
      end

      subroutine Hydrodyn_accel
      call Pressure_accel
      call Result_accel
      end

      subroutine Result_accel
      INCLUDE 'SPH_common.h'
      do j=1,dim
        do i_=1,n_p(timebin)
          i = p_list(i_,timebin)
          if (i_list(i).and.gas(i)) then
            accel(i,j) = q_acc(i,j) + P_acc(i,j)
            accel(i,j) = accel(i,j) + g(i,j)
          endif
        end do
      end do
      call Coriolis
      end
