#ifndef _N_BODY
/* Friday the 12th of July 1996 */
!///////////////////////////////////////////////////////////////////////////////
!///////////////////////   SPH Thermodynamics   ////////////////////////////////
!///////////////////////////////////////////////////////////////////////////////


#if (defined(_H_1)||defined(_H_2))&&!defined(_NON_CHEMICAL)
#  define _NON_CHEMICAL
#endif


#ifdef _PLUS_
#  define _FULL_CORRECTION_
#endif


#  ifdef _SPH_ISOTHERM
      subroutine Temp_
      INCLUDE 'SPH_common.h'
      parameter(Ro=8.154319e-3)     /* [e/m] K^-1 */
#    ifdef   _H_1
      freedom = 3
      w_m = 1/(XX+.25*YY)
#    elif defined (_H_2)
      freedom = 5
      w_m = 2/(XX+.5*YY)
#    endif
      R_=w_m/Ro * 2/freedom
      do i_=1,n_p(timebin)
        i = p_list(i_,timebin)
        Temp(i)= R_*u(i)
      end do
      end
#  endif _SPH_ISOTHERM




      subroutine Pressures
      INCLUDE 'SPH_common.h'
#  ifdef _GALAXY_2D
      parameter(Ro=8.702500e-5)     /* [e/m] K^-1 */
#  endif
#ifdef _VERBOSE
      print*,'Pressures'
#endif
#  ifdef _NON_CHEMICAL
#    if   defined(_H_1)
      freedom = 3
      w_m = 1/(XX+.25*YY)
#    elif defined (_H_2)
      freedom = 5
      w_m = 2/(XX+.5*YY)
#    endif
      call Monaghan_tensor
#  endif _NON_CHEMICAL
      do i_=1,n_p(timebin)
        i = p_list(i_,timebin)
        if(gas(i).and.i_list(i)) then
#  if !(defined(_GALAXY_2D)||defined(NON_CHEMICAL))
          !temperature --> T, P, 1/mu, HI, degrees of freedom
          Temp(i)=Temperature(u(i),rho(i),P(i),w_m,HI,freedom)
#  else
          Temp(i)=(2/freedom)*w_m*u(i)/Ro
          P(i) = (2/freedom)*u(i)/rho(i)
#  endif
#  ifdef _VERBOSE
c         print*,u(i),rho(i)*P(i),w_m,HI,freedom
c         print*,'T=',Temp(i),i
#  endif
        endif
      end do
#ifdef _VERBOSE
      print*,'Pressures: done'
#endif
      end


#ifndef _SPH_ISOTHERM
      subroutine SPH_first_law
      INCLUDE 'SPH_common.h'
#  if defined(_TUBE) || defined(_SPH_ADIABAT)
      call Pressures        !--> Temp, P
#  else
      call Cooling_rate !--> Temp, P, Lambda
#  endif
      call Monaghan_tensor
      call Viscosity_heating
      call Adiabatic_heating
#  if !defined(_SPH_ADIABAT)
      call Total_rate
#  endif
      end


#  if !( defined (_SPH_ADIABAT) )
      subroutine Total_rate
      INCLUDE 'SPH_common.h'
      do i_=1,n_p(timebin)
        i = p_list(i_,timebin)
        if(gas(i).and.i_list(i))u_dot(i) = lambda(i) + u_dot(i)
      end do
      end
#  endif


      subroutine adiabatic_heating
      INCLUDE 'SPH_common.h'
#ifdef _VERBOSE
      print*,'# adiabatic_heating'
#endif
      do i_=1,n_p(timebin)
        i = p_list(i_,timebin)
        if(gas(i).and.i_list(i))then
#  ifdef _FULL_CORRECTION_
          u_corr = 0
#  endif
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
#  ifdef _FULL_CORRECTION_
              u_corr=u_corr+
     &        0.25*mass(j)*(P(i)+P(j)+QQ(i,j_))*AddTerm(i,j)
#  endif
          end do
          u_left = u_left * P(i)
          u_dot(i) = .5*(u_left + u_rght)+u_dot_visc(i)
#  ifdef _FULL_CORRECTION_
     &                   + u_corr
#  endif
        endif
      end do
#ifdef _VERBOSE
      print*,'# adiabatic_heating: done'
#endif
      end


      subroutine viscosity_heating
      INCLUDE 'SPH_common.h'
#ifdef _VERBOSE
      print*,'# viscosity_heating'
#endif
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
#ifdef _VERBOSE
      print*,'# viscosity_heating: done'
#endif
      end

#endif


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
#define D0      5.1955e+4       /* Kelvin */
#define Ro      8.154319e-3     /* [e/m] K^-1 */
#define XX      0.75            /*  Typical     */
#define YY      0.25            /*   abundance  */
!     this procedure calculates the adiabatic speed of sound:

            ! gamma - 1 = ( P/rho )/u
            ! c^2 = gamma ( gamma - 1 ) u
            !     = (P/rho) [1 + P/(rho u)]

      INCLUDE 'SPH_common.h'
#ifdef _VERBOSE
      print*,'# sound'
#endif
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(gas(i).and.i_list(i))then
          p_by_rho=P(i)*rho(i)
#if !( defined (_H_2) || defined (_H_1) )
          /*
           Gambearra 1:
          */
          ! cs(i) = u(i)+p_by_rho
          /*
          Gambearra 2:
          */
          cs(i) = Ro*(3*(XX+.5*YY)*Temp(i)+D0*XX)*0.5
#else
          /*
            This is the usual to perform the adiabatic
            speed of sound:
          */
          !print*,i,u(i)
          cs(i) = p_by_rho * ( 1 + p_by_rho/u(i) )
#endif
          if (cs(i).lt.0) then
            print*,flsh,'error: invalid speed of sound!',norm
            print*,'c^2 =',cs(i)
            stop
          endif
#ifdef _SPMHD
          /*
          Gambearra 3:
          Magnetic case:
          c^2 -> c^2 + B^2/(8pi rho)
          -------------
          first method:
          -------------
          */
          /*
          B2=0
          do k=1,dim
            B2=B2+B(i,k)*B(i,k)
          end do
          */
          /*
          --------------
          second method:
          --------------
          */
          v2=0
          do k=1,dim
            v2=v2+v(i,k)*v(i,k)
          end do
          B2=0
          do k=1,dim
            Bk=0
            do j=1,dim
            do l=1,dim
              /* Bk is the perpendicular component of B onto v: */
              Bk=Bk+levi_civita(k,j,l)*B(i,j)*v(i,l)
            end do
            end do
            B2=B2+Bk*Bk
          end do
          B2=B2/(v2+0.0005)
          cs(i) = cs(i) + B2/(8*pi*rho(i)+0.0005)
#endif
          cs(i)=sqrt(cs(i))
        endif
      end do
#ifdef _VERBOSE
      print*,'# sound: done'
#endif
      end



      subroutine Monaghan_tensor
      INCLUDE 'SPH_common.h'
      call Sound_speed
#ifdef _VERBOSE
      print*,'# Monghan Tensor'
#endif
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
#ifdef _VERBOSE
      print*,'# Monaghan tensor: done.'
#endif
      end
#endif
