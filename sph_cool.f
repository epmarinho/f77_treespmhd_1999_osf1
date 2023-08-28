!///////////////////////////////////////////////////////////////////////////////
!////////////////   COOLING AND HEATING SECTION:   /////////////////////////////
!///////////////////////////////////////////////////////////////////////////////

#if !(defined(_SPH_ADIABAT)||defined(_N_BODY)||defined(_SPH_ISOTHERM)||defined(_GALAXY_2D))

      subroutine cooling_rate
      real Lamb,n_H,n_H1
      INCLUDE 'SPH_common.h'
        common/H_H2_CO/ Lamb(0:1023,0:63),Hy(0:1023,0:63)
      freedom = 5
      do il=1,n_p(timebin)
        i = p_list(il,timebin)
        if (gas(i).and.i_list(i)) then
          Temp(i) = Temperature(u(i),rho(i),P(i),w_m,HI,freedom)
#ifdef _VERBOSE
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
#endif
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
#if defined(_DAMP_COOLING_RATE)&&!(defined(_SPH_ADIABAT)||defined(_N_BODY)||defined(_SPH_ISOTHERM))
      if (you_can) call damp_cooling_rate
#endif
      end

#elif defined(_GALAXY_2D)&&!defined(_SPH_ADIABAT)

      subroutine cooling_rate
      INCLUDE 'SPH_common.h'
      u0=1.3e+0
      tau0=0.1 != 10^7 yr
      rate=1/tau0
      do il=1,n_p(timebin)
        i = p_list(il,timebin)
        if(gas(i).and.i_list(i))then
          LAMBDA(i) = -rate*(u(i)-u0)
        endif
      end do
      if (you_can) call damp_cooling_rate
      end

#endif

#if defined(_DAMP_COOLING_RATE)&&!(defined(_SPH_ADIABAT)||defined(_N_BODY)||defined(_SPH_ISOTHERM))

      subroutine damp_cooling_rate
      integer il,i
      real a,c
      INCLUDE 'SPH_common.h'
#ifdef _VERBOSE
      print*,'* damp_cooling_rate'
#endif
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
#ifdef _VERBOSE
      print*,'* damp_cooling_rate: done'
#endif
      end

#endif
