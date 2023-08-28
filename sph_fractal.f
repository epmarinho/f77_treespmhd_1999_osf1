      program SPH_fractal

        include 'SPH_common.h'
        !** print *,' file for initial conditions:'
        !** read(*,'(a132)')name
        name='BBmag_01_gas1.d'
        open(1,file=name,status='old')
        n=0
        i=0
    5     i=i+1
          read(1,*,end=10)
     &           mass(i)
     &          ,(xp(i,j),j=1,dim),(v(i,j),j=1,dim)
     &          ,u(i)
     &          ,(B(i,j),j=1,dim)
     &          ,rho(i)
        go to 5
   10   n=i-1
   11   close(1)
        first=.true.
        call maketree
        n_fix=48
        do i=1,n
          call smoothing_length(i)
        end do
        open(1,file='SPH_fractal.dat',status='unknown')
        do i=1,n
          write(1,'(8e16.8)')
     &      mass(i)
     &     ,(xp(i,j),j=1,dim)
     &     ,rho(i)
     &     ,h(i) 
     &     ,-log(rho(i))
     &     ,log(h(i))
        end do
        close(1)

      end
