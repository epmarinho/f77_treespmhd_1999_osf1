      PROGRAM sph_Tree_Gravitation
!
!
!      (individualy varying particle time-steps)
!
!
!
!
!
!
!
!
!
!
!      Author:      Eraldo Pereira Marinho
!
!
!
!
!
!
!
!                  Instituto Astronomico e Geofisico - USP.
!
!                  (011) 577 8599 (branch 235)
!
!                  eraldo@astro1.iagusp.usp.br
!
!                  eraldo@iag.usp.ansp.br
!
!                  FPSP::in%IAGUSP::ERALDO
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!      Version: 24.3 28-5-1994  
!
!
!            Last changings: 
!
!
!
!
!
!      Notes:
!
!            Physical unities are given in
!
!            [l] = 1 pc = 3.086 e+18 cm
!
!            [t] = 1 Myr = 3.156 e+13 s
!
!            [m] = 2.224 e+2 solar masses = 442.35 e+33 g
!            
!
!
!
!
!/////////////////////////// Main program //////////////////////////////
!
!23456789012345678901234567890123456789012345678901234567890123456789012
!        1         2         3         4         5         6         7
!
      print*,'////////////////////////////////////////////'
      print*
#ifdef _SPH_PLUS_
      print*,'Version SPH_PLUS'
#endif
      print*
      print*,'////////////////////////////////////////////'
      print*
      print*,'Vector/Parallel SPH treecode v. 11-12-1995'
      print*,'      Version 5.0'
      print*
      print*
      print*,'Designed by ERALDO PEREIRA MARINHO'
      print*
      print*,'      Instituto Astronomico e Geofisico.'
      print*,'      Universidade de Sao Paulo - Brazil'
      print*
      print*,'(P) 1992. All Rights Reserved.'
      print*
      print*,'////////////////////////////////////////////'
      print*
      print*
      call dosimulation
      print *
      print *,'Done!'
!
      end
!///////////////////////////////////////////////////////////////////////

