#ifndef _N_BODY
!///////////////////////////////////////////////////////////////////////////////
!////////////////////////   SPH USEFUL SUBROUTINES   ///////////////////////////
!///////////////////////////////////////////////////////////////////////////////
/*
    Last revisions:

       Thursday the 7th of the July 1996

       Sunday the 4th of the August 1996
*/


      subroutine Densities
      /*
         Densities are calculated by the combined gather-scatter technique
         due to Hernquist & Katz 1989. For that, it crucial the employement
         of a gather NNS followed by a scatter NNS.
      */
      INCLUDE 'SPH_common.h'
#ifdef _VERBOSE
      print*,'|-> densities'
#endif
      call Neighboring
      do i_=1,n_p(timebin)
        if(gas(p_list(i_,timebin)))call Rho_(p_list(i_,timebin))
      end do
#ifdef _VERBOSE
      print*,'|-> densities: done'
#endif
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
        /* smoothing kernel in gather-scatter form (HK89): */
        w_ij = .5*(sphKernel(r_ij,h(i))+sphKernel(r_ij,h(j)))
        rho(i) = w_ij * mass(j) + rho(i)
      end do
#if !( defined (_TUBE) || defined (_SPH_ADIABAT) )
      /*
      For realm ISM simulations only:
      impose a density threshold of 1.2e-4 (~1 cm^-3)
      */
      if(rho(i).lt.(1.2e-4))rho(i)=1.2e-4
#else
      if(rho(i).eq.0)then
            print*,'error: zero density'
            print*,'      particle',i
            print*,'      ',neighb(i),' neighbors'
            print*,'      h =',h(i)
            print*,'      time-bin',timebin
            stop
      endif
#endif
      end



      subroutine Get_kernel_grad
      INCLUDE 'SPH_common.h'
      common /aux/ aux(nn,dim)
      equivalence (aux,x_ij)
      dimension x_ij(nn,dim)
      /* for timebin particles do: */
#ifdef _VERBOSE
      print*,'kernel gradients'
#endif
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
            /*
               it is assumed the gather-scatter approach 
               grad_ww is the polar derivative:
                  dW_ij/dr_ij
               with r_ij being the separation between particles i and j
            */
              grad_ww=.5*(D_kernel(r_ij,h(i))+D_kernel(r_ij,h(j)))
            /*
               now, it is calculated the vector gradient with respect
               to i-particle's position
            */
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
#ifdef _VERBOSE
      print*,'kernel gradients: done!'
#endif _VERBOSE
#  ifdef _SPH_PLUS_
      !print*,'Calculating kernel gradients corrections...'
      call Get_Modif_grad
#  endif
      end



#ifdef _NEW_VERSION
      function V_divergence(i)
        INCLUDE 'SPH_common.h'
        /* extract rho*div(v) for h(i) */
        h_dot(i) = 0
        do j_=1,neighb(i)
          j=neighb_list(i,j_)
          if (gas(j)) then
            do k=1,dim
              h_dot(i)=h_dot(i)+(v(j,k)-v(i,k))*mass(j)*grad_w(i,j_,k)
            end do
          endif
        end do
        /* convert rho*div(v) to h_dot */
        h_dot(i) = h_dot(i) * tau / rho(i)
      end
#endif


      subroutine Get_h_dot
      /*
      Warning:
            Check method before starting a serious simulation !
      */
      INCLUDE 'SPH_common.h'
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        if(h_list(i).and.gas(i))then
          /* extract rho*div(v) for h(i) */
          h_dot(i) = 0
          do j_=1,neighb(i)
            j=neighb_list(i,j_)
            if (gas(j)) then
              do k=1,dim
                h_dot(i)=h_dot(i)+(v(j,k)-v(i,k))*mass(j)*grad_w(i,j_,k)
              end do
            endif
          end do
          /* mean half time-step */
          tau=.25*(dt(timebin)+dt_old(i))
          /* convert rho*div(v) to h_dot */
          h_dot(i) = h_dot(i) * tau / (rho(i) * 3)
        endif
      end do
      end




#ifdef _CHECKUP
!///////////////////////////////////////////////////////////////////////////////
!//////////////////////////   SECURITY CHECKINGS   /////////////////////////////
!///////////////////////////////////////////////////////////////////////////////



      subroutine Discard_unknown(i)
      INCLUDE 'SPH_common.h'
      do itime=n_timebins,1,-1
      i_=1
      do while(i_.le.n_p(itime))
        ip=p_list(i_,itime)
        if (ip.eq.i) then
          !print*,'particle',i,' discarded'
          j=i_
          do while(j.lt.n_p(itime))
            p_list(j,itime)=p_list(j+1,itime)
            j=j+1
          end do
          n_p(itime)=n_p(itime)-1
          n=n-1
          !print*,'remaining ',n
          return
        else
          i_=i_+1
        endif    
      end do
      end do
      end


      subroutine Discard_neighb(i,k)
      INCLUDE 'SPH_common.h'
      i_=1
      do while (i_.le.neighb(k))
      ip=neighb_list(k,i_)
      if (ip.eq.i) then
        j=i_
        do while (j.lt.neighb(k))
          neighb_list(k,j)=neighb_list(k,j+1)
          j=j+1
        end do
        neighb(k)=neighb(k)-1
        return
      else
        i_=i_+1
      endif
      end do
      end


      subroutine Discard(i)
      INCLUDE 'SPH_common.h'
      i_=1
      !print*,n_p(timebin)
      do while (i_.le.n_p(timebin))
        ip=p_list(i_,timebin)
        if(ip.eq.i)then
          !print*,'particle',i,' discarded'
          j=i_
          do while(j.lt.n_p(timebin))
            p_list(j,timebin)=p_list(j+1,timebin)
            j=j+1
          end do
          n_p(timebin)=n_p(timebin)-1
          n=n-1
          !print*,'remaining ',n
          return
        else
          i_=i_+1
        endif    
      end do
      end
#endif
#endif
