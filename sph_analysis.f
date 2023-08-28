      Program SPH_analysis
        print*,'Performing SPH_analysis'
        call Perform_Analysis
        print*,' READY'
      end




      subroutine Assembl_SPH
        character*14 err_distrib,analysis,sphfile
        common/default_files/ err_distrib,analysis
        common/divergence2/ div_min,div_max,r_min,r_max
        INCLUDE 'SPH_common.h'
        common/div/     	div_r(nn)
        common/normalization/     y_b(nn,dim)
        print*,'Assembl_SPH:'
        print*
00010        print*,'enter option:'
        print*,'    (1): read .core file'
        print*,'    (2): read .d file'
        print*,'    (3): generate a cubic latice'
        print*,'    (4): quit'
        read(*,*)iop
        if(iop.eq.1) then
          err_distrib='other_dist.dat'
          analysis   ='other_ana.dat'
          call readcore
        else
        if(iop.eq.2) then
          err_distrib='other_dist.dat'
          analysis   ='other_ana.dat'
          print *,' enter file name:'
          read(*,'(a132)')name
          !name=name(1:index(name,' ')-1)
          write(*,*)'the file name is',name,index(name,' ')
          print *,' enter N (0 for the maximum of the file, noting that'
     &    ,' the default',nn,' is the maximum):'
          read(*,*)n
          if((n.eq.0).or.(n.gt.nn))n=nn
          call read_phase
        else
        if(iop.eq.3) then
          err_distrib='latic_dist.dat'
          analysis   ='latic_ana.dat'
          call assembl_particles
        else
        if(iop.eq.4)then
          stop
        else
          goto 10
        endif
        endif
        endif
        endif
        sphfile='sphdata.dat'
        /****************************************
        start it up:
        *****************************************/
        timebin=1
        n_timebins=1
        n_p(timebin)=n
        print*,'n =',n
        do i=1,n_p(timebin)
          gas(i)=.true.
          i_list(i)=gas(i)
          p_list(i,timebin)=i
        end do
#ifndef _TUBE
        n_fix = 32
#else
        n_fix = 4
#endif
        first=.true.
        hydrosim = .true.
        print*,'Performing essential calculations:'
        print*,'Maketree:'
        call maketree
        tol_nfix = 1
        print*,'Search_lengths:'
        call search_lengths
        print*,'Densities:'
        call densities
        do i=1,n
          if (.not.(rho(i).gt.0))then
            print*,'Error: zero density. Check Desnsities'
            stop
          endif
        end do
        open(1,file=sphfile,status='unknown')
        do i=1,n
          write(1,'(2e15.7)') rho(i), h(i)
        end do
        close(1)
        print*,'Get_kernel_grad:'
        call get_kernel_grad
#ifndef _SPMHD
        print*,'Calc_h_dot:'
        call Calc_h_dot
        print*,'Mass_transport:'
        call mass_transport
        print*,'done.'
        print*,'it was created a file with the analysis of mass conservation;'
        print*,'wish you to continue with the rest of calculations?'
        read(*,'(a1)')answer
        if(answer.eq.'n')stop
        call OtherCalcs
#endif
      end


      subroutine OtherCalcs
        INCLUDE 'SPH_common.h'
#ifndef _SPMHD
        common/normalization/     y_b(nn,dim)
#endif
        common/divergence2/ div_min,div_max,r_min,r_max
        common/div/     	div_r(nn)
        r_min=0
        r_max=2

#ifndef _SPMHD
        /* Clear error gradient */
        print*,'# 01'
        do i=1,n
          do j=1,dim
            y_b(i,j)=0
          end do
        end do
#endif
#ifdef _SPMHD
        div_min=0
        div_max=div_min
#else
        div_min=0
        div_max=6
#endif
        print*,'calculating normalization, gradient and divergence'
        do i=1,n
          div_r(i) = 0
#ifndef _SPMHD
#endif
          do jnb=1,neighb(i)
            j=neighb_list(i,jnb)
            /* write(*,'(a4,i5)'),'j = ',j */
            r_ij=0
            do k=1,dim
#ifdef _SPMHD
              B_ij =B(i,k)- B(i,k)
              div_r(i) = div_r(i) + B_ij * mass(j) * grad_w(i,jnb,k)
#else
              x_ji=xp(j,k)-xp(i,k)
              r_ji=r_ji+x_ji*x_ji
              div_r(i) = div_r(i) + x_ji * mass(j) * grad_w(i,jnb,k)
              y_b(i,k) = y_b(i,k) + mass(j) * grad_w(i,jnb,k)/rho(j)
#endif
            end do
#ifndef _SPMHD
            w_ij = .5*(sphKernel(r_ij,h(i))+sphKernel(r_ij,h(j)))
            r_ij=sqrt(r_ij)
#endif
          end do
          div_r(i) = div_r(i) / rho(i)
#ifndef _SPMHD
#endif

#ifdef _SPMHD
          div_min=amin(div_min, div_r(i))
          div_max=amax(div_max, div_r(i))
#endif

/*
          do jnb=1,neighb(i)
            j=neighb_list(i,jnb)
            write(*,'(a4,i5)'),'j = ',j
            do k=1,dim
              x_ij=xp(j,k)/(rho(j)*rho(j))+xp(i,k)/(rho(i)*rho(i))
              div_r(i) = div_r(i) + x_ij * mass(j) * grad_w(i,jnb,k)
            end do
          end do
          div_r(i) = div_r(i) * rho(i)
*/
        end do
        print*,'done.'

      end



      Subroutine Perform_Analysis
        character*14 err_distrib,analysis
        common/default_files/ err_distrib,analysis
        parameter (ngrd=48)
        common/divergence2/     div_min,div_max,r_min,r_max
        INCLUDE 'SPH_common.h'
        common/div/    	div_r(nn)
        common/list/    	list(nn)
        common/normalization/     y_b(nn,dim)
        common/histog/ f(0:ngrd)
        ! logical cutoff

        call Assembl_SPH
#ifdef _CUTOFF
        print*,'wish you to cut off edge particles?'
        read(*,'(a)')answer
        cutoff=(answer.eq.'y')
        if(cutoff)then
          r_min=.85
          div_min=2.5
        endif
#endif

        scale=ngrd/(div_max-div_min)
        scale2=ngrd/(r_max-r_min)

        do i=0,ngrd
          f(i)  = 0
        end do

        open(1,file=analysis,status='unknown')
        print*,'scale 1:',scale,' scale 2:',scale2
        do i=1,n
          igrd=int(scale*(div_r(i)-div_min))
          !print*,igrd
          if((igrd.ge.0).and.(igrd.le.ngrd))then
            f(igrd)=f(igrd)+1
          endif
#ifndef _SPMHD
          write(1,'(9e15.7)')
     &              div_r(i)		! 1
     &            , rho(i)		! 2
     &            ,(y_b(i,j),j=1,dim)	! 3 - 5
     &            ,(xp(i,j), j=1,dim)	! 6 - 8
     &            , h(i)			! 9
#else
          write(1,'(2e15.7)')
     &              div_r(i)		! 1
     &            , rho(i)		! 2
     
#endif
        end do
        close(1)
        tot_f=0
        tot_g=0
        print*,'normalizing for ',ngrd,' grid points'
        do igrd=0,ngrd
          tot_f=tot_f+f(igrd)
          f(igrd)=-f(igrd)*scale
        end do
        print*,'Total f=',tot_f
c        print*,'Total g=',tot_g
        tot_f=-tot_f
#ifndef _SPMHD
c        tot_g=-tot_g
#endif
        tot = 0
        do igrd=0,ngrd
          f (igrd) = f (igrd)/tot_f
          tot = tot + f(igrd)
        end do
        !print*,'tot_f =',tot_f
        print*,'integral check:',tot/scale,' done.'

        ! smoothing ou by median
        ! print*,'smoothing out result'
        ! do igrd=0,ngrd-1
          ! f(igrd)=.5*(f(igrd)+f(igrd+1))
        ! end do
        ! print*,'done.'

        ! output result
        open(1,file=err_distrib,status='unknown')
        do igrd=0,ngrd
#ifndef _SPMHD
          write(1,*)igrd/scale+div_min,f(igrd)
#else
          write(1,*)igrd/scale+div_min,f(igrd)
#endif
        end do
        close(1)

        return
      end



      Subroutine assembl_particles
        external ran
        INCLUDE 'SPH_common.h'
        isemente=1236527631
        jsemente=2001235445
        ksemente=1423423143
        print*,'Assembl_particles:'
        print*,'enter latice parameter:'
        read(*,*)rlatice
        print*,'enter number of linear grid points:'
        read(*,*)l
        gl = (float(nn)**(1./3) - 1)
        l  = l - 1
        if (l.lt.1) then
            l = int (gl)
        else
            l  = int ( amin ( gl, float(l) ) )
        endif
        print*,'# of linear grid points:',l+1
        n=0
        do i=0,l
          do j=0,l
            do k=0,l
              n=n+1
              mass(n)=1.0
              xp(n,1) = float (k)/l + (rlatice*(.5-ran(ksemente)) )
              xp(n,2) = float (j)/l + (rlatice*(.5-ran(jsemente)) )
              xp(n,3) = float (i)/l + (rlatice*(.5-ran(isemente)) )
            end do
          end do
        end do
        return
      end



      Subroutine Calc_h_dot
        INCLUDE 'SPH_common.h'
            							      /*
          Calculates the smoothed particle velocity-divergence multiplied by
          rho:

        print*,'Calc_h_dot:'
        print*
            							      */
        do i=1,n
          h_dot(i) = 0
          do jnb=1,neighb(i)
            j=neighb_list(i,jnb)
            do k=1,dim
              v_ij = v(i,k)-v(j,k)
              h_dot(i) = h_dot(i)  +  v_ij*mass(j)*grad_w(i,jnb,k)
            end do
          end do
              							      /*
          h_dot has the physical dimension of the velocity-divergence
          but it is divided by 3 since it is a linear scale:

            h_dot/h = 1/3 div (v)
              							      */
          if (rho(i).eq.0)then
            print*,'rho null for',i
            stop
          else
            h_dot(i) = - (h_dot(i) / (rho(i) * 3))
          endif
        end do
        return
      end



      function W_h(i,j,l)
        INCLUDE 'SPH_common.h'
        r_ij = 0
        do k=1,dim
          x_ij = xp(i,k)-xp(j,k)
          r_ij = r_ij +  x_ij*x_ij
        end do
        r_ij = sqrt(r_ij)
            							      /*
          W_h is the W derivative with respect to h, multiplied by h.
          The following expression works only for K(r)/h^3 kernels [v.
          Marinho (1996)].
            							      */
        W_h = -(r_ij*D_kernel(r_ij,h(l))+3*sphKernel(r_ij,h(l)))
        return
      end



      Subroutine mass_transport
        INCLUDE 'SPH_common.h'
        print*,'Mass_transport:'
        open(1,file='mass_transp.dat',status='unknown')
        do i=1,n
          varpi = 0
          do jnb=1,neighb(i)
            j=neighb_list(i,jnb)
            do k=1,dim
              varpi = varpi +
     &            0.5*mass(j)*(
     &              h_dot(i)*W_h(i,j,i) +
     &              h_dot(j)*W_h(i,j,j)
     &            )
            end do
          end do
          write(1,'(3e15.7,i4)')
     &          rho(i),varpi,3*h_dot(i)*rho(i),neighb(i)
        end do
        close(1)
      end
