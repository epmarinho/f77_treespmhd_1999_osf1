# 1 "SPH_integral.F"
! 
! 
! 


      subroutine Integral_quantities
      INCLUDE 'SPH_common.h'
      call Clear_integrals
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
      INCLUDE 'SPH_common.h'
      ! clear kinetic-energy accumulator:
      T_tot = 0
      do j=1,dim
        do i=1,n
          v_j = .5*((2-factor)*v(i,j)+factor*v_o(i,j))
          v2 = v_j * v_j
          T_tot = T_tot + .5*v2*mass(i)
        end do
      end do
      end

