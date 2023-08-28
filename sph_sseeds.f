#ifndef _N_BODY
      subroutine StartSeeds
!
      INCLUDE 'SPH_common.h'
!
      open(5,file='SPH_seeds.dat',status='old')
      read(5,*)iroot
      close(5)
      do i=1,n
            iseed(i) = 1791123513*ran(iroot)
            iseed(i)=iseed_max(iseed(i))
            seed=seed + iseed(i)/float(n)
            jseed(i) = 2008235129*ran(iroot)
            jseed(i)=iseed_max(jseed(i))
            seed=seed + jseed(i)/float(n)
            kseed(i) = 1235213459*ran(iroot)
            kseed(i)=iseed_max(kseed(i))
            seed=seed + kseed(i)/float(n)
            lseed(i) = 1611235341*ran(iroot)
            lseed(i)=iseed_max(lseed(i))
            seed=seed + lseed(i)/float(n)
            mseed(i) = 1352631411*ran(iroot)
            mseed(i)=iseed_max(mseed(i))
            seed=seed + mseed(i)/float(n)
            nseed(i) = 1251336413*ran(iroot)
            nseed(i)=iseed_max(nseed(i))
            seed=seed + nseed(i)/float(n)
      end do
      iroot=int(seed/6.)
      iroot=iseed_max(iroot)
      open(5,file='SPH_seeds.dat',status='unknown')
      !write(5,*)iroot
      close(5)
      return
      end

      function iseed_max(seed)
      integer seed
      seed=max(seed,2008235129-seed)
      seed=abs(float(seed))
      seed=seed+1-mod(seed,2)
      iseed_max=seed
      return
      end
#endif

