# 1 "SPH_intschm.F"





      subroutine predictor_corrector (dE)
      INCLUDE 'SPH_common.h'
      if(n_p(timebin).lt.1)return
      print*,'predictor_corrector'
      call initial_saving
      call move_particles
      call set_no_convergence
      call reset_time_list
      call predict_velocities(1.,0)
      call Preliminaries
      call extern_gravity_2d
      call integrate_on_u_v
      call save_old_timesteps
      print*,'predictor_corrector: done'
      end

      subroutine Integrate_on_u_v
      logical repeat
      repeat= .true.
      print*,'Integrate_on_u_v'
      do while (repeat)
        call Step_counter
        call predict_velocities(1.,1)
        call advance_energies
        call advance_velocities
        call check_convergence(repeat)
      end do
      print*,'Integrate_on_u_v: done'
      end

      subroutine Set_i_list
      INCLUDE 'SPH_common.h'
      do i_=1,n_p(timebin)
        i_list(p_list(i_,timebin))=.true.
      end do
      end

      subroutine Check_convergence(repeat)
      logical repeat
      INCLUDE 'SPH_common.h'
      itol=23
      tol=1
      do i=1,itol
        tol=.5*tol
      end do
      repeat = .false.
      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        if(gas(i).and.i_list(i))then
          if(u(i).gt.0)then
            control=abs(u(i)-u_old(i))-tol*u(i)
            i_list(i)=(control.gt.0)
            control=(abs(u(i)-u_old(i))-tol*u_old(i))
            i_list(i)=i_list(i).or.(control.gt.0)
          endif
          do j=1,dim
            control=abs(v(i,j)-v_old(i,j))-tol*v(i,j)
            i_list(i)=i_list(i).or.(control.gt.0)
            control=abs(v(i,j)-v_old(i,j))-tol*v_old(i,j)
            i_list(i)=i_list(i).or.(control.gt.0)
          end do
          repeat = (repeat.or.i_list(i)).and.(iter.lt.128)
        endif
      end do
      end

      subroutine Preliminaries
      INCLUDE 'SPH_common.h'
      !print*,'call Maketree'
      call Maketree
      !print*,'call Search_lengths'
      call Search_lengths
      !print*,'call Densities'
      call Densities
      !print*,'call Get_kernel_grad'
      call Get_kernel_grad
      end

      subroutine Step_counter
      INCLUDE 'SPH_common.h'
      iter = iter + 1
      end

      subroutine Set_no_convergence
      integer il,i
      INCLUDE 'SPH_common.h'
      iter = 0
      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        i_list(i)=.true.
      end do
      end

      subroutine Initial_saving
      integer il,i,j
      INCLUDE 'SPH_common.h'
      do j=1,dim
        do il=1,n_p(timebin)
          i=p_list(il,timebin)
          v_o(i,j)=v(i,j)
        end do
      end do
      do il=1,n_p(timebin)
        i=p_list(il,timebin)
        if(gas(i))then
          tau=.25*(dt(timebin)+dt_old(i))
          u_o(i)=u(i)+tau*u_dot(i)! u_dot is the old value
        endif
      end do
      end

      subroutine Save_old_velocities
      integer il,i
      INCLUDE 'SPH_common.h'
      do il=1,n_p(timebin)
       i=p_list(il,timebin)
       do j=1,dim
         if(gas(i).or.i_list(i)) v_old(i,j) = v(i,j)
       end do
      end do
      end

      subroutine Save_old_energies
      integer il,i
      INCLUDE 'SPH_common.h'
      do il=1,n_p(timebin)
       i=p_list(il,timebin)
       if(gas(i).or.i_list(i)) u_old(i) = u(i)
      end do

      end




      subroutine initialize
      INCLUDE 'SPH_common.h'
      if(not_yet)then
	print*,'initialize: start'



        n_p(1)      = n
        n_timebins  = 1
        timebin     = 1
        do i=1,n_p(timebin)
          p_list(i,timebin)=i

          i_list(i) = .true.
          gas(i)    = .true.
        end do
        call Preliminaries
        call Sound_speed
        call Pressures



        call Get_h_dot
        you_can  = (ichoice.eq.0)




        starting = .false.
        not_yet  = .false.
        first    = .false.
        if(dt(0).eq.0.0)return
      else
        starting = .true.
      endif
      if(ichoice.eq.1)then

        call SPH_First_law


        call Accelerations







        call Start_list
        do i=1,n
          dt_old(i) = 0
        end do
        ichoice=0
        you_can = .true.
        call predictor_corrector(0)
      endif
      call Start_list

      call Set_i_list

      print*,'initialize: done'
      end



      subroutine Integration_scheme (dE)
      INCLUDE 'SPH_common.h'
      print*,'Integration_scheme'

      call StartSeeds

      call Initialize
      if(dt(0).eq.0.0)return
      call multi_leapfrog (dE)



      call integral_quantities
      call recover_dtm
      print*,'Integration_scheme:done'
      end
