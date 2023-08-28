	character*256 nameout
        logical isotherm
        integer dim,dn
	real mass
!
!       parameter(nn=8192,dim=3)
!       parameter(nn=16384,dim=3)
!       parameter(nn=24576,dim=3)
!	parameter(nn=32768,dim=3)
	parameter(nn=262144,dim=3)

!
        common /matrices/ x(nn,dim),v(nn,dim)
        common /vectors/ u(nn),r(nn)
	common /field/ B(dim)
        common /scalars/ gamma,omega,pi,rho_m
        common /integers/ iseed,jseed,kseed,n
        common /booleans/ isotherm
	common /strings/ nameout
