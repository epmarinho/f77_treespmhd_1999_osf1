      subroutine vorticity
        include 'SPH_common.h'
	dimension omega(nn,dim)
	logical ciclica

        do i=1,dim
          do j=1,dim
            do k=1,dim
              if((i*j*k).eq.6)then
                l=i+3*(j+3*k)
                ciclica=(l.eq.20).or.(l.eq.24).or.(l.eq.34)
                if(ciclica)then
                  levi_civita(i,j,k) =  1
                else
                  levi_civita(i,j,k) = -1
                endif
              else
                levi_civita(i,j,k)   =  0
              endif
            end do
          end do
        end do

        print*,'data-file name (with B-field):'
#ifdef _DEBUG
        name='BBmag_11_ortho_gas64.d'
#else
	read(*,'(a128)'),name
#endif
	open(1,file=name,status='old')
        i=0
    5     i=i+1
          read(1,*,end=10)
     &      mass(i)
     &     ,(xp(i,j),j=1,dim),(v(i,j),j=1,dim)
     &     ,u(i)
     &     ,(B(i,j),j=1,dim)
     &     ,rho(i)
          if(i.eq.n)go to 11
        go to 5
   10     n=i-1
   11   close(1)
        print *
        print *,n,' phase-points were read'
        print *

	first=.true.
	call maketree

	print*,'calculating smoothing lengths'
	do i=1,n
C	  call Particle_length (i)
	  call smoothing_length (i)
	end do

	print*,'searching for neighborings'
        do i=1,n
          neighb(i) = 0
          call set_permit
          do node=1,maxnode
            if(permit(node))then
              call effective_neighbors(i,2*h(i))
            else
              call forbid(node)
            endif
          end do
        end do
C	do i=1,n
C          neighb(i) = 0
C          call Set_permit
C          do node=1,maxnode
C            if(permit(node))then
C              call SRN(i,2*h(i))
C            else
C              call Forbid(node)
C            endif
C          end do
C	end do

	print*,'calculating kernel gradient'
	do i=1,n
          do j_=1,neighb(i)
            j=neighb_list(i,j_)
            r_ij = 0
            do k=1,dim
              r_ij = r_ij + (xp(i,k)-xp(j,k))*(xp(i,k)-xp(j,k))
            end do
            r_ij=sqrt(r_ij)
            if(r_ij.gt.0)then
              grad_ww=.5*(D_kernel(r_ij,h(i))+D_kernel(r_ij,h(j)))
              do k=1,dim
                grad_w(i,j_,k) =(xp(i,k)-xp(j,k))*grad_ww / r_ij
              end do
            else
              do k=1,dim
                grad_w(i,j_,k) = 0
              end do
            endif
          end do
	end do

	print*,'calculating vorticity'
	do i=1,n
	  do k=1,dim
	    omega(i,k)=0
	  end do
	  do j_=1,neighb(i)
	    j=neighb_list(i,j_)
	    do k=1,dim
	      cross=0
	      do i1=1,dim
	        do i2=1,dim
	          cross=cross+
     &	          grad_w(i,j_,i1)*(v(j,i2)-v(i,i2))*levi_civita(i1,i2,k)
                end do
	      end do
	      omega(i,k)=omega(i,k)+
     &	      mass(j)*cross
            end do
	  end do! j_
	end do! i
	do i=1,n
	  do k=1,dim
	    omega(i,k)=omega(i,k)/rho(i)
	  end do
	end do! i

	print*,'outputing data'
	open(1,file='vorticity.d',status='unknown')
	write(1,'(10e16.8)')(
     &   mass(i)		!1
     &  ,(xp(i,j),j=1,dim)	!2
     &  ,(v(i,j),j=1,dim)	!5
     &  ,(omega(i,j),j=1,dim)	!8
     &  ,i=1,n)
        close(1)

      end
      
      program calc_vorticity
        call vorticity
      end
