# 1 "SPH_usr.F"


      subroutine Densities
      INCLUDE 'SPH_common.h'
      print*,'|-> densities'
      call Neighboring
      do i_=1,n_p(timebin)
        if(gas(p_list(i_,timebin)))call Rho_(p_list(i_,timebin))
      end do
      print*,'|-> densities: done'
      end

      subroutine Rho_(i)
      INCLUDE 'SPH_common.h'
      rho(i)=0
      do j_=1,neighb(i)
        j=neighb_list(i,j_)
        r_ij = 0
        do k=1,dim
          x_ij = xp(i,k)-xp(j,k)
          r_ij = r_ij + x_ij*x_ij
        end do
        r_ij = sqrt(r_ij)
        w_ij = .5*(sphKernel(r_ij,h(i))+sphKernel(r_ij,h(j)))
        rho(i) = w_ij * mass(j) + rho(i)
      end do
      if(rho(i).lt.(1.2e-4))rho(i)=1.2e-4
      end

      subroutine Get_kernel_grad
      INCLUDE 'SPH_common.h'
      common /aux/ aux(nn,dim)
      equivalence (aux,x_ij)
      dimension x_ij(nn,dim)
      do i_=1,n_p(timebin)
        i = p_list(i_,timebin)
        if(gas(i))then
          do j_=1,neighb(i)
            j=neighb_list(i,j_)
            r_ij = 0
            do k=1,dim
              x_ij(1,k) = xp(i,k) - xp(j,k)
              r_ij = r_ij + x_ij(1,k)*x_ij(1,k)
            end do
            r_ij=sqrt(r_ij)
            if(r_ij.gt.0)then
              grad_ww=.5*(D_kernel(r_ij,h(i))+D_kernel(r_ij,h(j)))
              do k=1,dim
                grad_w(i,j_,k) = x_ij(1,k) * grad_ww / r_ij
              end do
            else
              do k=1,dim
                grad_w(i,j_,k) = 0
              end do
            endif
          end do
        endif
      end do
      end

      subroutine Get_h_dot
      INCLUDE 'SPH_common.h'
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(h_list(i).and.gas(i))then
          h_dot(i) = 0
          do j_=1,neighb(i)
            j=neighb_list(i,j_)
            if (gas(j)) then
              do k=1,dim
                h_dot(i)=h_dot(i)+(v(j,k)-v(i,k))*mass(j)*grad_w(i,j_,k)
              end do
            endif
          end do
          tau=.25*(dt(timebin)+dt_old(i))
          h_dot(i) = h_dot(i) * tau / (rho(i) * 3)
        endif
      end do
      end
