      subroutine test_parallelization(n,dim,y,r)
        parameter(epsilon=0.01)
        integer n,i,j,dim
        real x(n,dim),y(n,dim),r
        do i=1,n
          r=epsilon
          do j=1,dim
            r=r+x(i,j)*x(i,j)
          end do
          r=sqrt(r)
          do j=1,dim
            y(i,j)=-x(i,j)/(r*r*r)
          end do
        end do
      end
