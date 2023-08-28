      subroutine extern_gravity_2D
        include 'SPH_common.h'
        real Mcore

        Mcore=1e+0
        Rcore=1e-1

        do il=1,n_p(timebin)
          i=p_list(il,timebin)
          r=Rcore*Rcore
          do j=1,dim
            r=r+xp(i,j)*xp(i,j)
          end do
          r=sqrt(r)
          factor=-Mcore/(r*r*r)
          !factor=0
          do j=1,dim
            g(i,j)=xp(i,j)*factor
          end do
          !*** write(*,'(4e16.7)')(xp(i,j),j=1,dim),(g(i,j),j=1,dim)
        end do

        !*** stop

      end
