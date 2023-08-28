!///////////////////////////////////////////////////////////////////////////////
!/////////////////////   INTEGRAL QUANTITIES   /////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////
      subroutine Integral_quantities
      INCLUDE 'SPH_common.h'
#if defined (_TUBE) || defined (_NON_SELFGRAVITATING)
      E_tot=0
#else
      call Treepotential
      E_tot = W_tot
#endif
#if !defined (_TUBE)
      call Angular_momentum
#endif
      call Kinetic_energy(1.)
      E_tot = E_tot + T_tot
#if !defined (_N_BODY)
      call Thermal_energy
      E_tot = E_tot + U_tot
      call Pressure_virial
#  if defined (_SPMHD)
      call Magnetic_energy
      E_tot = E_tot + E_mag
#  endif
#endif
      end

#ifndef _TUBE
      subroutine Clear_Ang_Momentum
      INCLUDE 'SPH_common.h'
      do j=1,dim
        angmom(j)=0
      end do
      end

      subroutine Angular_momentum
      INCLUDE 'SPH_common.h'
      call Clear_Ang_Momentum
      do j=1,dim
        vmoment(j) = 0
!        for all system particles do:
        do i=1,n
#ifndef _N_BODY
          v_j = .5*(v(i,j)+v_o(i,j))
#else
          v_j = v(i,j)-.5*g(i,j)*dt_old(i)
#endif
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
#endif

#ifndef _N_BODY
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
#ifdef _VERBOSE
      print*,'Thermal_energy =',U_tot
#endif
      end

!Check whether or not magnetic stress contribution is missing!!!
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
#endif

      subroutine Kinetic_energy(factor)
      !Don't ask me why I have placed the magnetic energy calculations here!
      INCLUDE 'SPH_common.h'
      ! clear kinetic-energy accumulator:
      T_tot = 0
      do j=1,dim
        do i=1,n
#if defined (_N_BODY)
          v_j = v(i,j)-.5*g(i,j)*dt_old(i)
#else
          v_j = .5*((2-factor)*v(i,j)+factor*v_o(i,j))
#endif
          v2 = v_j * v_j
          T_tot = T_tot + .5*v2*mass(i)
        end do
      end do
      end

#if defined (_SPMHD)
      subroutine Magnetic_energy
      E_mag = 0
      do j=1,dim
        do i=1,n
#   ifdef _BACKGROUND_FIELD_
          Bij_2=B(i,j)+Bbg(j)
#   else
          Bij_2=B(i,j)
#   endif
          Bij_2=Bij_2*Bij_2
          E_mag=E_mag+Bij_2*mass(i)/rho(i)
        end do
      end do
      E_mag = E_mag/(8*pi)
      end
#endif
