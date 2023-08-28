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
      B_rate_k=0
      do jnb=1,neighb(i)
        j=neighb_list(i,jnb)
        ! Bscal = B^ dot grad_w^
        Bscal=0
        vscal = 0
c       performs the scalar product (v(j,l)-v(i,l)) * grad_w(i,jnb,l)
c       and the scalar product B(i,l)*grad_w(i,jnb,l)
        do l=1,dim
          vscal = vscal + (v(j,l)-v(i,l)) * grad_w(i,jnb,l) !Ok!
          Bscal = Bscal + B(i,l)*grad_w(i,jnb,l) !Ok!
        end do
        div      = div      + vscal*mass(j) !Ok!
        B_rate_k = B_rate_k+(Bscal*mass(j))*(v(j,k)-v(i,k)) !Ok!
      end do
      B_rate_k=(B_rate_k-B(i,k)*div)/(rho(i)+0.0000005)
c     return result:
      B_rate = B_rate_k
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
