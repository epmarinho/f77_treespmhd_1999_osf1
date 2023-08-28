# 1 "SPH_mhd.F"



! Function BB_grad(i,j,k) takes the k-component of the SPH-gradient from the
! Maxwell's stress-tensor.

      function BB_grad(i,jnb,j)
      INCLUDE 'SPH_common.h'
      B2=0
      do k=1,dim
        B2=B2+B(i,k)*B(i,k)
      end do
      Bscal=0
      do k=1,dim
        Bscal=Bscal+B(i,k)*grad_w(i,jnb,k)
      end do
      BB_grad=(B2*grad_w(i,jnb,j)-2*B(i,j)*Bscal)/
     &(8*pi*rho(i)*rho(i))
      end




      subroutine field_contrib
      INCLUDE 'SPH_common.h'
      do k=1,dim
        do i_=1,n_p(timebin)
          i=p_list(i_,timebin)
          if(gas(i))then
            XB(i,k)=0
            do jnb=1,neighb(i)
              j=neighb_list(i,jnb)
              XB(i,k)=XB(i,k) + mass(j)*(BB_grad(i,jnb,k)
     &                        + BB_grad(j,jnb,k))
            end do
            XB(i,k) = -XB(i,k)
          endif
        end do
      end do
      call background_field
      end


      subroutine background_field
      include 'SPH_common.h'
      !print*,'background_field'
      do k=1,dim
        do i_=1,n_p(timebin)
          i=p_list(i_,timebin)
          if(gas(i))then
            contrib=0
            contrib1=0
            contrib2=0
            contrib3=0
            contrib4=0
            !scl2 -> Bbg^ dot B^_i/rho_i^2
            scl2=0
            do l=1,dim
              scl2=scl2+Bbg(l)*B(i,l)
            end do
            scl2=scl2/(rho(i)*rho(i))
            !for each neighbor do:
            do jnb=1,neighb(i)
              j=neighb_list(i,jnb)
              !common scalar products:
              !scl1 -> Bbg^ dot Grad^_i W_ij
              !scl3 -> Bbg^ dot B^_j
              scl1=0
              scl3=0
              do l=1,dim
                scl1=scl1+Bbg(l)*grad_w(i,jnb,l)
                scl3=scl3+Bbg(l)*B(j,l)
              end do! l
              !contrib1 -> Sum_j m_jBbg^dotGrad^_iW_ij B_j/rho_j^2(k-direction)
              contrib1=contrib1+
     &                  mass(j)*scl1*
     &                  B(j,k)/(rho(j)*rho(j))
              !contrib2 -> Sum_j m_j Bbg^dot Grad^_iW_ij
              contrib2=contrib2+mass(j)*scl1
              !contrib3 -> Sum_j m_j Grad^_iW_ij Bbg^dot B_j/rho_j^2
              contrib3=contrib3+
     &                  mass(j)*grad_w(i,jnb,k)*
     &                  scl3/(rho(j)*rho(j))
              !contrib4 -> Sum_j m_j Grad^_iW_ij (k-direction)
              contrib4=contrib4+
     &                  mass(j)*grad_w(i,jnb,k)
            end do! each neighbor
            contrib2=contrib2*B(i,k)/(rho(i)*rho(i))
            contrib4=contrib4*scl2
            contrib=(contrib1+contrib2-(contrib3+contrib4))/(4*pi)
            XB(i,k)=XB(i,k)+contrib
          endif
        end do
      end do
      !print*,'background_field: end'
      end
