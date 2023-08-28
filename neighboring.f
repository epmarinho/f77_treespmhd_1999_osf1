      subroutine neighboring
      include 'SPH_common.h'
      do i_=1,n_p(timebin)
        i=p_list(i_,timebin)
        neighb(i) = 0
        call set_permit
        do node=1,maxnode
          if(permit(node))then
            call effective_neighbors(i,2*h(i))
          else
            call forbid
          endif
        end do
      end do
      end
