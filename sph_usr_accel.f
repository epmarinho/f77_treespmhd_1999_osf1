#ifndef _N_BODY
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
#  ifdef _BACKGROUND_PRESSURE_
#    ifndef _GALAXY_2D
              P_k        = P_k + gw_ijk*(P(j)-Pbg/(rho(j)*rho(j)))
#    else
              P_k        = P_k + gw_ijk*(P(j)-Pbg/rho(j))
#    endif
#  else
              P_k        = P_k + gw_ijk*P(j)
#  endif
              P_acc(i,k) = P_acc(i,k) + gw_ijk
            end do
            q_acc(i,k) = -q_acc(i,k)
#  ifdef _BACKGROUND_PRESSURE_
#    ifndef _GALAXY_2D
            P_acc(i,k) = -(P_k+P_acc(i,k)*(P(i)-Pbg/(rho(i)*rho(i))))
#    else
            P_acc(i,k) = -(P_k+P_acc(i,k)*(P(i)-Pbg/rho(i)))
#    endif
#  else
            P_acc(i,k) = -(P_k+P_acc(i,k)*P(i))
#  endif
          endif
        end do
      end do
#ifdef _TUBE
      call Boundary
      end

      subroutine Boundary
      INCLUDE 'SPH_common.h'
      do k=1,dim
        do i=1,neighb(1)+2
          P_acc(i,k)=0
        end do
        do i=n-(2+neighb(n)),n
          P_acc(i,k)=0
        end do
      end do
      do i=1,neighb(1)+2
          P(i)=P(neighb(1))
      end do
      do i=n-(2+neighb(n)),n
          P(i)=P(n-neighb(n))
      end do
#endif
      end


      subroutine Accelerations /* Called from Preliminars only */
#ifndef _TUBE
#  ifndef _NON_SELFGRAVITATING
      call treegravity
#  else
      call extern_gravity_2d
#  endif
#endif
      call Hydrodyn_accel
      end



      subroutine Hydrodyn_accel
      call Pressure_accel
#ifdef _SPMHD
      call Field_contrib
#endif
      call Result_accel
      end

 

      subroutine Result_accel
      INCLUDE 'SPH_common.h'
      do j=1,dim
        do i_=1,n_p(timebin)
          i = p_list(i_,timebin)
          if (i_list(i).and.gas(i)) then
                   accel(i,j) = q_acc(i,j) + P_acc(i,j)
#ifndef _TUBE
#ifdef _SPMHD
     &               + XB(i,j)
#endif
            accel(i,j) = accel(i,j) + g(i,j)
#endif
          endif
        end do
      end do
#     ifdef _ROTATING_FRAME
      call Coriolis
#     endif _ROTATING_FRAME
      end
#endif
