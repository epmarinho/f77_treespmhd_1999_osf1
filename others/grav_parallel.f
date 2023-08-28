C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 19:09:07
# 1 "SPH_grav.F"
 
! 
!                  TREE-GRAVITY SUBROUTINES
! 
 
      SUBROUTINE TREEGRAVITY
       INTEGER CELL, IL, I, J, L
!     INCLUDE 'SPH_common.h'
# 1 "SPH_common.F"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      SYSTEM DEFAULTS FOR THE SPH-TREECODE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      Author:  ERALDO PEREIRA MARINHO.
!
!                  Instituto Astronomico e Geofisico - USP.
!
!                  (011) 577 8599
!
!                  eraldo@astro1.iagusp.usp.br
!
!                  eraldo@iag.usp.ansp.br
!
!                  IAGUSP::ERALDO
!
!      Version: 4.0 10-31-1995
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!      data-type declaration:
!
!
 
       INTEGER DIM, ONE, TWO, THREE, FOUR, EIGHT, DOWN, T, COUNT, ROW, 
     X   NODE, NTBINS, TIME_LEVEL, GLOBAL_TIME_LEVEL, TIMEBIN, P_LIST, 
     X   LABEL, VE, DEEPEST, NN, MAXOCTANTS, NNODE, MAXNEIGHB, ICHOICE, 
     X   IS, NS, LST, IFILE, NI, KWSC, N_TIMEBINS, NTBINS2, N, NEXT, 
     X   MAXNODE, N_P, NP, IFIRST, IISEED, KRONECKER, LEVI_CIVITA, 
     X   STTMNT_COUNT, ITER, ITERCOUNTER, ITERPLUS, ITER_NEIGHB, N_FIX, 
     X   NEIGHB_LIST, NEIGHB, ISEED, JSEED, KSEED, LSEED, MSEED, NSEED, 
     X   TOL_NFIX
 
 
 
 
 
       REAL MASS, CELLMASS, EPSILON, G, AREST, PI, EPSILON2, DT_ELAPSED,
     X    DT, D, DT_OLD, P_INERTIA, X_CELL, ANGMOM, S, V, XP, THETA, 
     X   THETA2, TOT_MASS, E_TOT, W_TOT, T_TOT, U_TOT, E_MIN, E_MAX, EO,
     X    COURANT, VMOMEN, LAMBDA, LUMINOSITY, TWO_THIRD, XNORM, ROOT_3,
     X    SC, W_P, W_S, CC, ACCEL, RHO, V_O, V_OLD, U, U_DOT, U_O, H, 
     X   H_DOT, T1, T2, UCC, U_MASS, U_LENG, FINE, U_TIME, U_DENS, 
     X   U_DEN_N, XX, YY, ZZ, DTM_SAFE, ALPHA, BETA, ETA, GUESSING, P, 
     X   QQ, P_ACC, Q_ACC, GRAD_W, U_DOT_VISC, TEMP, CS, U_OLD, VMON
 
 
 
 
 
 
 
 
 
       LOGICAL NEXTNODE, PERMIT, CONVERGENCE, YOU_CAN, STARTING, 
     X   CONVERGING, OK, NONZERO, NONOVER, NZ, NV, OVC, NOT_YET, FIRST, 
     X   H_LIST, I_LIST, GAS
 
       CHARACTER ESC, CR, CUTL*4, REV*4, FLSH*4, HLTD*4, NORM*4, HOME*4,
     X    CLS*4, ANSWER
       CHARACTER*132 EVOLUTION, NAME, NAMEOUT, NAMECORE, RELAT, PERCENT
 
 
 
 
 
       REAL*8 U_ENER, U_VOL
 
!
!      global parameters:
!
!
       PARAMETER (PI = 3.141592654, TWO_THIRD = .6666667)
 
 
 
 
 
       PARAMETER (XNORM = .3183098861)
 
       PARAMETER (ROOT_3 = 1.732050808)
       PARAMETER (CC = 0.744438)
       PARAMETER (FINE = 0.56)
!
!            stars parameters:
 
 
 
 
 
 
 
       PARAMETER (T1 = 3.37, T2 = 4.34, UCC = 1.343295)
!
!
!      system array sizes:
!      -------------------
 
 
 
 
 
 
 
       PARAMETER (DIM = 3, VERT = 8.)
 
 
 
 
 
 
 
 
 
 
 
       PARAMETER (NN = 40960, MAXOCTANTS = 327672, NNODE = 40959)
 
 
 
 
 
 
 
 
 
 
!
!      usual integers and maximum number of neighbors:
!      -----------------------------------------------
       PARAMETER (ZERO = 0., ONE = 1, TWO = 2, THREE = 3, FOUR = 4, 
     X   EIGHT = 8)
 
 
 
       PARAMETER (MAXNEIGHB = 96)
 
 
 
 
 
 
!      time resolution:
!      ----------------
 
 
 
 
 
 
 
       PARAMETER (NTBINS = 4)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
CGS --> computer units:
 
 
 
 
       PARAMETER (U_MASS = 0.2248201E-35, U_LENG = 0.3138732E-18, U_TIME
     X    = 0.3168568E-13, U_DENS = 0.7270638E+20, U_DEN_N = 
     X   0.1207747E-03, U_ENER = DBLE (0.2206062E-45), U_VOL = DBLE (
     X   0.3092165E-55))
 
Cloud abundances:
       PARAMETER (XX = .75, YY = .25, ZZ = 0.)
 
 
 
 
!      COMMON VARIABLES
 
!            PART I: STRINGS
 
!      File names:
 
       COMMON /ANSWER/ ANSWER
       COMMON /CHR/ ESC, CR, CUTL, REV, FLSH, HLTD, NORM, HOME, CLS
       COMMON /NAMES/ EVOLUTION, NAME, NAMEOUT, NAMECORE, RELAT, PERCENT
 
 
 
!
!            PART II: 32-BITS SCALARS
 
 
 
 
!
       COMMON /CHOICE/ ICHOICE
       COMMON /IFNEWFLAG/ FIRST, YOU_CAN, NOT_YET
       COMMON /RELATE/ IS, NS, IFILE, NI, INI, NIT, ITER, ITERCOUNTER, 
     X   ITERPLUS, ITER_NEIGHB
       COMMON /NEXTNODE/ NEXTNODE
       COMMON /CONVERGE/ CONVERGENCE, CONVERGING, OK, NONZERO, NONOVER, 
     X   NZ, NV, OVC
       COMMON /STARTING_NOW/ STARTING
       COMMON /WELL_SEP_NODE_COUNTER/ KWSC
       COMMON /TIME_LEVELS/ TIME_LEVEL, GLOBAL_TIME_LEVEL
 
       COMMON /TIME_BINS/ TIMEBIN, N_TIMEBINS, NTBINS2, DEEPEST, DTM, 
     X   DTM_SAFE
       COMMON /FIXED_NEIGHBOR_NUMBER/ N_FIX
       COMMON /NEIGHB_GUESSING/ GUESSING
       COMMON /NEIGHB_TOLER/ SC, TOL_NFIX
       COMMON /NEIGHB_Y/ Y
       COMMON /EXTRAS1/ ALPHA, BETA, ETA
 
 
 
 
 
       COMMON /EXTRAS2/ THETA, EPSILON, THETA2, EPSILON2
 
       COMMON /NODE/ NODE, MAXNODE
       COMMON /NEXT/ NEXT
       COMMON /TOTAL_NUMBER/ N
       COMMON /TOTAL_MASS/ TOT_MASS
 
 
 
 
       COMMON /ENERGIES/ E_TOT, E_MIN, E_MAX, EO, W_TOT, T_TOT, U_TOT, 
     X   W_S, W_P, LUMINOSITY
 
       COMMON /COURANT_FACTOR/ COURANT
       DIMENSION ROW(2)
 
 
 
 
!
!
!
!            PART III: SINGLE 3-D VECTORS
!
!
       COMMON /ROWS/ ROW
       DIMENSION ANGMOM(3), VMOMENT(3)
       COMMON /MOMENTA/ VMOMENT, RR3, ANGMOM
       DIMENSION KRONECKER(3,3)
 
 
 
       COMMON /IDENTITY/ KRONECKER
       DIMENSION LEVI_CIVITA(3,3,3)
       COMMON /LEVI/ LEVI_CIVITA
       DIMENSION IFIRST(8)
       COMMON /IFIRST/ IFIRST
       DIMENSION COUNT(8)
       COMMON /COUNT/ COUNT
       DIMENSION H_LIST(40960)
!
!
!            PART IV: Rn-VECTORS
!
 
       COMMON /HLIST/ H_LIST
       DIMENSION GAS(40960)
       COMMON /STATUS/ GAS
       DIMENSION MASS(40960)
 
       COMMON /PARTICLE_MASS/ MASS
       DIMENSION S(40960)
       COMMON /PARTICLE_SIZE/ S
       DIMENSION NEIGHB(40960)
 
       COMMON /NUMBER_OF_NEIGHBOURS/ NEIGHB
       DIMENSION H_DOT(40960), H(40960)
       COMMON /SMOOTHING_LENGTHS/ H, H_DOT
       DIMENSION RHO(40960)
       COMMON /SPH_DENSITIES/ RHO
       DIMENSION P(40960)
       COMMON /PRESSURE/ P
       DIMENSION TEMP(40960)
       COMMON /TEMPERATURES/ TEMP
       DIMENSION CS(40960)
       COMMON /SPEED_OF_THE_SOUND/ CS
       DIMENSION VMONAG(40960), LAMBDA(40960), U_DOT_VISC(40960), U_DOT(
     X   40960), U_OLD(40960), U_O(40960), U(40960)
       COMMON /THERMAL_ENERG/ U, U_O, RR4, U_OLD, U_DOT, RR5, U_DOT_VISC
     X   , LAMBDA, RR6, VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, RR7, KSEED, LSEED, RR8, MSEED, NSEED
       DIMENSION I_LIST(40960)
       COMMON /NOT_CONV_LIST/ I_LIST
 
       COMMON /TREE_SEED/ IISEED
       DIMENSION T(40960,4)
!
!       data configuration:
!          nn=n
!            nnode = n - 1
!            maxoctants = vert * nnode
!
!        the array t has the following aspect:
!      -------------------------------------
!                    /t/
!      -------------------------------------
!        /  row  /|           j
!      -------------------------------------
!            i  |  1  |  2  |   3  |   4   
!      -------------------------------------
!            1  |  t  |label| cell | iflg  
!            2  |  "  | "   |  "   |   "   
!            3  |  "  | "   |  "   |   "   
!             ... | ... |...  | ...  |  ...  
!              n  |  "  | "   |  "   |   "   
!
!      where, t = vert * t + cell is a peano parameter, which 
!      distinguish a node from other:
!             cell is the octant number in {1,2,...,vert};
!             iflg = 1 if particle is a leaf, otherwise it equals 0;
!             label is the fixed particle index.
!
       COMMON /T/ T
       DIMENSION NEIGHB_LIST(40960,96)
!
!
!            PART V: Rmn-VECTORS
!
 
       COMMON /SPH_LISTS/ NEIGHB_LIST
       DIMENSION QQ(40960,96)
       COMMON /VISCOTENSOR/ QQ
       DIMENSION GRAD_W(40960,96,3)
 
 
 
 
       COMMON /GRAD_LIST/ GRAD_W
       DIMENSION P_ACC(40960,3), Q_ACC(40960,3), G(40960,3)
 
 
!
!
!            PART VI: Rn3 VECTORS
!
 
 
       COMMON /FORCES/ G, RR9, Q_ACC, RR10, P_ACC
       DIMENSION N_P(4)
 
 
 
 
 
 
 
 
 
 
!
!
!            PART VII: Rn nt VECTORS
!
!
       COMMON /TIME_BIN/ N_P
       DIMENSION P_LIST(40960,4)
       COMMON /TIMELIST/ P_LIST
       DIMENSION DT_OLD(40960)
       COMMON /INDIVIDUAL_TIME/ DT_OLD
       DIMENSION DT(0:4)
       COMMON /TIME_STEPS/ DT
       DIMENSION DT_ELAPSED(4)
       COMMON /ELAPSED_TIME/ DT_ELAPSED
       DIMENSION LST(40960,8)
 
 
 
!
!            PART VIII: TREE-STRUCTURES  
!
       COMMON /OCTLIST/ LST
       DIMENSION PERMIT(40959)
       COMMON /PERMIT/ PERMIT
       DIMENSION AREST(40959)
       COMMON /TREE_REAL_0/ AREST
       DIMENSION CELLMASS(40959,8)
       COMMON /TREE_REAL_1/ CELLMASS
       DIMENSION X_CELL(40959,8,3)
       COMMON /TREE_REAL_3/ X_CELL
       DIMENSION P_INERTIA(40959,8,3,3)
       COMMON /TREE_REAL_6/ P_INERTIA
       DIMENSION DOWN(40959,8), LABEL(40959,8), NP(40959,8)
       COMMON /TREE_INTG/ NP, RR11, LABEL, RR12, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, RR13, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
C*$* Padding ( RR13, RR12, RR11, RR10, RR9, RR8, RR7, RR6, RR5, RR4, RR3
C*$*&  )
       INTEGER II13, II12, II11, II10, II9, II8, II7
       INTEGER*2 HI2, HI1
       PARAMETER (II13 = 2, II12 = 4, II11 = 3, II10 = 2730, II9 = 1)
       PARAMETER (HI2 = 0, II8 = 8, II7 = 0, HI1 = 1)
       EXTERNAL MPPTID, MPTEPA, MPTXPA, MPPFKD, MPOFRK, PKTREEGRAVITY_
       INTEGER MPPTID, MPPFKD, II19(0:9)
       CHARACTER AA6*28
       INTEGER II18, II6, II5, II4, II3, II2, II1, NODE2, NODE1
       REAL RR2, RR1
       REAL RR13 (86), RR12 (86), RR11 (86), RR10 (86), RR9 (86), RR8 (
     X   86), RR7 (86), RR6 (86), RR5 (86), RR4 (86), RR3 (86)
       AUTOMATIC II18
       DATA II19(0)/0/ 
       DATA II19(1)/0/ 
       DATA II19(2)/28/ 
       DATA AA6/';SPH_grav.F;TREEGRAVITY;13;;'/ 
       EQUIVALENCE (AA6, II19(3))
       INTEGER II21, II20
       PARAMETER (II21 = 0, II20 = 500)
       II18 = MPPTID ()
      !dimension xa(dim),ga(dim)
       PRINT *, 'treegravity'
       DO IL=HI1,N_P(TIMEBIN)
        I = P_LIST(IL,TIMEBIN)
        IF (MAXNODE .GT. II20 .AND. MPPFKD () .EQ. II21) THEN
         CALL MPOFRK (PKTREEGRAVITY_,II19)
        ELSE
C!!!!! PARALLEL REGION IF (MAXNODE .GT. 500) SHARED (MAXNODE,PERMIT) 
C!!!!!& LOCAL (NODE1)
         CALL MPTEPA (II18)
C!!!!! PARALLEL DO 
         DO NODE1=HI1,MAXNODE
          PERMIT(NODE1) = .TRUE.
         END DO
C!!!!! END PARALLEL DO NOWAIT
         CALL MPTXPA (II18)
C!!!!! END PARALLEL REGION 
        END IF
        IF (MAXNODE .GT. II7) NODE = MAXNODE
        KWSC = II7
        ! Performing tree-walk
        II1 = MAXNODE
        DO NODE2=HI1,II1
         NODE = NODE2
          ! vectorization is broken here:
         IF (PERMIT(NODE2)) THEN
          DO CELL=HI1,II8
           IF (NP(NODE,CELL) .GT. HI2) CALL GRAVITY (I,CELL)
          END DO
         ELSE
          CALL FORBID
         END IF
        END DO
        NODE = MAX0 (II1, II7) + HI1
        II5 = MOD (KWSC - II9, II10) + II9
        !-----------------------------
        !do j=1,dim
        !  xa(j)=xp(i,j)
        !end do
        !call phi_galax(xa,ga,epsilon)
        ! Now, the particle has a list
        ! of kwsc cell'contributions:
        !-----------------------------
        DO J=HI1,II11
         G(I,J) = 0.
        END DO
        II3 = II9
        II2 = II5
        DO II6=II9,KWSC,II10
         II4 = II3 + II2 - II9
        !-----------------------------
        !do j=1,dim
        !  xa(j)=xp(i,j)
        !end do
        !call phi_galax(xa,ga,epsilon)
        ! Now, the particle has a list
        ! of kwsc cell'contributions:
        !-----------------------------
         DO J=HI1,II11
          RR1 = G(I,J)
          DO L=II3,II4-II11,II12
           RR1 = RR1 + GNODE(L,J)
           RR1 = RR1 + GNODE(L+II9,J)
           RR1 = RR1 + GNODE(L+II13,J)
           RR1 = RR1 + GNODE(L+II11,J)
            !g(i,j)=g(i,j)+ga(j)
          END DO
          G(I,J) = RR1
          RR2 = G(I,J)
          DO L=L,II4,II9
           RR2 = RR2 + GNODE(L,J)
          END DO
          G(I,J) = RR2
         END DO
         II3 = II4 + II9
         II2 = II10
        END DO
        !-----------------------------
       END DO
       PRINT *, 'done!'
       RETURN 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 19:09:07
      SUBROUTINE PKTREEGRAVITY_ ( MPPID, MPPNPR )
       AUTOMATIC II14
       INTEGER II14
       INTEGER MPPID
       INTEGER NNODE1
       PARAMETER (NNODE1 = 40959)
       INTEGER MPPNPR
       INTEGER MAXNODE1
       INTEGER*2 HI11
       PARAMETER (HI11 = 1)
       INTEGER NODE3
       AUTOMATIC II17
       INTEGER II17
       AUTOMATIC NODE4
       INTEGER NODE4
       AUTOMATIC II16
       INTEGER II16
       AUTOMATIC II15
       INTEGER II15
       LOGICAL PERMIT1(40959)
       COMMON /NODE/ NODE3, MAXNODE1
       COMMON /PERMIT/ PERMIT1
C!!!!! PARALLEL REGION IF (MAXNODE .GT. 500) SHARED (MAXNODE,PERMIT) 
C!!!!!& LOCAL (NODE1)
       INTEGER II22
       PARAMETER (II22 = 1)
       II16 = MAXNODE1 - HI11 + II22
       II17 = (II16 + MPPNPR - II22) / MPPNPR
       II14 = HI11 + MPPID * II17
       II15 = MIN (MAXNODE1, II14 + (II17 - II22))
C!!!!! PARALLEL DO 
       DO NODE4=II14,II15
        PERMIT1(NODE4) = .TRUE.
       END DO
C!!!!! END PARALLEL DO NOWAIT
C!!!!! END PARALLEL REGION 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 19:09:07
 
 
      SUBROUTINE GRAVITY ( I, CELL )
       INTEGER CELL
!     INCLUDE 'SPH_common.h'
# 1 "SPH_common.F"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      SYSTEM DEFAULTS FOR THE SPH-TREECODE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      Author:  ERALDO PEREIRA MARINHO.
!
!                  Instituto Astronomico e Geofisico - USP.
!
!                  (011) 577 8599
!
!                  eraldo@astro1.iagusp.usp.br
!
!                  eraldo@iag.usp.ansp.br
!
!                  IAGUSP::ERALDO
!
!      Version: 4.0 10-31-1995
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!      data-type declaration:
!
!
 
       INTEGER DIM, ONE, TWO, THREE, FOUR, EIGHT, DOWN, T, COUNT, ROW, 
     X   NODE, NTBINS, TIME_LEVEL, GLOBAL_TIME_LEVEL, TIMEBIN, P_LIST, 
     X   LABEL, VE, DEEPEST, NN, MAXOCTANTS, NNODE, MAXNEIGHB, ICHOICE, 
     X   IS, NS, LST, IFILE, NI, KWSC, N_TIMEBINS, NTBINS2, N, NEXT, 
     X   MAXNODE, N_P, NP, IFIRST, IISEED, KRONECKER, LEVI_CIVITA, 
     X   STTMNT_COUNT, ITER, ITERCOUNTER, ITERPLUS, ITER_NEIGHB, N_FIX, 
     X   NEIGHB_LIST, NEIGHB, ISEED, JSEED, KSEED, LSEED, MSEED, NSEED, 
     X   TOL_NFIX
 
 
 
 
 
       REAL MASS, CELLMASS, EPSILON, G, AREST, PI, EPSILON2, DT_ELAPSED,
     X    DT, D, DT_OLD, P_INERTIA, X_CELL, ANGMOM, S, V, XP, THETA, 
     X   THETA2, TOT_MASS, E_TOT, W_TOT, T_TOT, U_TOT, E_MIN, E_MAX, EO,
     X    COURANT, VMOMEN, LAMBDA, LUMINOSITY, TWO_THIRD, XNORM, ROOT_3,
     X    SC, W_P, W_S, CC, ACCEL, RHO, V_O, V_OLD, U, U_DOT, U_O, H, 
     X   H_DOT, T1, T2, UCC, U_MASS, U_LENG, FINE, U_TIME, U_DENS, 
     X   U_DEN_N, XX, YY, ZZ, DTM_SAFE, ALPHA, BETA, ETA, GUESSING, P, 
     X   QQ, P_ACC, Q_ACC, GRAD_W, U_DOT_VISC, TEMP, CS, U_OLD, VMON
 
 
 
 
 
 
 
 
 
       LOGICAL NEXTNODE, PERMIT, CONVERGENCE, YOU_CAN, STARTING, 
     X   CONVERGING, OK, NONZERO, NONOVER, NZ, NV, OVC, NOT_YET, FIRST, 
     X   H_LIST, I_LIST, GAS
 
       CHARACTER ESC, CR, CUTL*4, REV*4, FLSH*4, HLTD*4, NORM*4, HOME*4,
     X    CLS*4, ANSWER
       CHARACTER*132 EVOLUTION, NAME, NAMEOUT, NAMECORE, RELAT, PERCENT
 
 
 
 
 
       REAL*8 U_ENER, U_VOL
 
!
!      global parameters:
!
!
       PARAMETER (PI = 3.141592654, TWO_THIRD = .6666667)
 
 
 
 
 
       PARAMETER (XNORM = .3183098861)
 
       PARAMETER (ROOT_3 = 1.732050808)
       PARAMETER (CC = 0.744438)
       PARAMETER (FINE = 0.56)
!
!            stars parameters:
 
 
 
 
 
 
 
       PARAMETER (T1 = 3.37, T2 = 4.34, UCC = 1.343295)
!
!
!      system array sizes:
!      -------------------
 
 
 
 
 
 
 
       PARAMETER (DIM = 3, VERT = 8.)
 
 
 
 
 
 
 
 
 
 
 
       PARAMETER (NN = 40960, MAXOCTANTS = 327672, NNODE = 40959)
 
 
 
 
 
 
 
 
 
 
!
!      usual integers and maximum number of neighbors:
!      -----------------------------------------------
       PARAMETER (ZERO = 0., ONE = 1, TWO = 2, THREE = 3, FOUR = 4, 
     X   EIGHT = 8)
 
 
 
       PARAMETER (MAXNEIGHB = 96)
 
 
 
 
 
 
!      time resolution:
!      ----------------
 
 
 
 
 
 
 
       PARAMETER (NTBINS = 4)
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
CGS --> computer units:
 
 
 
 
       PARAMETER (U_MASS = 0.2248201E-35, U_LENG = 0.3138732E-18, U_TIME
     X    = 0.3168568E-13, U_DENS = 0.7270638E+20, U_DEN_N = 
     X   0.1207747E-03, U_ENER = DBLE (0.2206062E-45), U_VOL = DBLE (
     X   0.3092165E-55))
 
Cloud abundances:
       PARAMETER (XX = .75, YY = .25, ZZ = 0.)
 
 
 
 
!      COMMON VARIABLES
 
!            PART I: STRINGS
 
!      File names:
 
       COMMON /ANSWER/ ANSWER
       COMMON /CHR/ ESC, CR, CUTL, REV, FLSH, HLTD, NORM, HOME, CLS
       COMMON /NAMES/ EVOLUTION, NAME, NAMEOUT, NAMECORE, RELAT, PERCENT
 
 
 
!
!            PART II: 32-BITS SCALARS
 
 
 
 
!
       COMMON /CHOICE/ ICHOICE
       COMMON /IFNEWFLAG/ FIRST, YOU_CAN, NOT_YET
       COMMON /RELATE/ IS, NS, IFILE, NI, INI, NIT, ITER, ITERCOUNTER, 
     X   ITERPLUS, ITER_NEIGHB
       COMMON /NEXTNODE/ NEXTNODE
       COMMON /CONVERGE/ CONVERGENCE, CONVERGING, OK, NONZERO, NONOVER, 
     X   NZ, NV, OVC
       COMMON /STARTING_NOW/ STARTING
       COMMON /WELL_SEP_NODE_COUNTER/ KWSC
       COMMON /TIME_LEVELS/ TIME_LEVEL, GLOBAL_TIME_LEVEL
 
       COMMON /TIME_BINS/ TIMEBIN, N_TIMEBINS, NTBINS2, DEEPEST, DTM, 
     X   DTM_SAFE
       COMMON /FIXED_NEIGHBOR_NUMBER/ N_FIX
       COMMON /NEIGHB_GUESSING/ GUESSING
       COMMON /NEIGHB_TOLER/ SC, TOL_NFIX
       COMMON /NEIGHB_Y/ Y
       COMMON /EXTRAS1/ ALPHA, BETA, ETA
 
 
 
 
 
       COMMON /EXTRAS2/ THETA, EPSILON, THETA2, EPSILON2
 
       COMMON /NODE/ NODE, MAXNODE
       COMMON /NEXT/ NEXT
       COMMON /TOTAL_NUMBER/ N
       COMMON /TOTAL_MASS/ TOT_MASS
 
 
 
 
       COMMON /ENERGIES/ E_TOT, E_MIN, E_MAX, EO, W_TOT, T_TOT, U_TOT, 
     X   W_S, W_P, LUMINOSITY
 
       COMMON /COURANT_FACTOR/ COURANT
       DIMENSION ROW(2)
 
 
 
 
!
!
!
!            PART III: SINGLE 3-D VECTORS
!
!
       COMMON /ROWS/ ROW
       DIMENSION ANGMOM(3), VMOMENT(3)
       COMMON /MOMENTA/ VMOMENT, RR7, ANGMOM
       DIMENSION KRONECKER(3,3)
 
 
 
       COMMON /IDENTITY/ KRONECKER
       DIMENSION LEVI_CIVITA(3,3,3)
       COMMON /LEVI/ LEVI_CIVITA
       DIMENSION IFIRST(8)
       COMMON /IFIRST/ IFIRST
       DIMENSION COUNT(8)
       COMMON /COUNT/ COUNT
       DIMENSION H_LIST(40960)
!
!
!            PART IV: Rn-VECTORS
!
 
       COMMON /HLIST/ H_LIST
       DIMENSION GAS(40960)
       COMMON /STATUS/ GAS
       DIMENSION MASS(40960)
 
       COMMON /PARTICLE_MASS/ MASS
       DIMENSION S(40960)
       COMMON /PARTICLE_SIZE/ S
       DIMENSION NEIGHB(40960)
 
       COMMON /NUMBER_OF_NEIGHBOURS/ NEIGHB
       DIMENSION H_DOT(40960), H(40960)
       COMMON /SMOOTHING_LENGTHS/ H, H_DOT
       DIMENSION RHO(40960)
       COMMON /SPH_DENSITIES/ RHO
       DIMENSION P(40960)
       COMMON /PRESSURE/ P
       DIMENSION TEMP(40960)
       COMMON /TEMPERATURES/ TEMP
       DIMENSION CS(40960)
       COMMON /SPEED_OF_THE_SOUND/ CS
       DIMENSION VMONAG(40960), LAMBDA(40960), U_DOT_VISC(40960), U_DOT(
     X   40960), U_OLD(40960), U_O(40960), U(40960)
       COMMON /THERMAL_ENERG/ U, U_O, RR8, U_OLD, U_DOT, RR9, U_DOT_VISC
     X   , LAMBDA, RR10, VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, RR11, KSEED, LSEED, RR12, MSEED, 
     X   NSEED
       DIMENSION I_LIST(40960)
       COMMON /NOT_CONV_LIST/ I_LIST
 
       COMMON /TREE_SEED/ IISEED
       DIMENSION T(40960,4)
!
!       data configuration:
!          nn=n
!            nnode = n - 1
!            maxoctants = vert * nnode
!
!        the array t has the following aspect:
!      -------------------------------------
!                    /t/
!      -------------------------------------
!        /  row  /|           j
!      -------------------------------------
!            i  |  1  |  2  |   3  |   4   
!      -------------------------------------
!            1  |  t  |label| cell | iflg  
!            2  |  "  | "   |  "   |   "   
!            3  |  "  | "   |  "   |   "   
!             ... | ... |...  | ...  |  ...  
!              n  |  "  | "   |  "   |   "   
!
!      where, t = vert * t + cell is a peano parameter, which 
!      distinguish a node from other:
!             cell is the octant number in {1,2,...,vert};
!             iflg = 1 if particle is a leaf, otherwise it equals 0;
!             label is the fixed particle index.
!
       COMMON /T/ T
       DIMENSION NEIGHB_LIST(40960,96)
!
!
!            PART V: Rmn-VECTORS
!
 
       COMMON /SPH_LISTS/ NEIGHB_LIST
       DIMENSION QQ(40960,96)
       COMMON /VISCOTENSOR/ QQ
       DIMENSION GRAD_W(40960,96,3)
 
 
 
 
       COMMON /GRAD_LIST/ GRAD_W
       DIMENSION P_ACC(40960,3), Q_ACC(40960,3), G(40960,3)
 
 
!
!
!            PART VI: Rn3 VECTORS
!
 
 
       COMMON /FORCES/ G, RR13, Q_ACC, RR14, P_ACC
       DIMENSION N_P(4)
 
 
 
 
 
 
 
 
 
 
!
!
!            PART VII: Rn nt VECTORS
!
!
       COMMON /TIME_BIN/ N_P
       DIMENSION P_LIST(40960,4)
       COMMON /TIMELIST/ P_LIST
       DIMENSION DT_OLD(40960)
       COMMON /INDIVIDUAL_TIME/ DT_OLD
       DIMENSION DT(0:4)
       COMMON /TIME_STEPS/ DT
       DIMENSION DT_ELAPSED(4)
       COMMON /ELAPSED_TIME/ DT_ELAPSED
       DIMENSION LST(40960,8)
 
 
 
!
!            PART VIII: TREE-STRUCTURES  
!
       COMMON /OCTLIST/ LST
       DIMENSION PERMIT(40959)
       COMMON /PERMIT/ PERMIT
       DIMENSION AREST(40959)
       COMMON /TREE_REAL_0/ AREST
       DIMENSION CELLMASS(40959,8)
       COMMON /TREE_REAL_1/ CELLMASS
       DIMENSION X_CELL(40959,8,3)
       COMMON /TREE_REAL_3/ X_CELL
       DIMENSION P_INERTIA(40959,8,3,3)
       COMMON /TREE_REAL_6/ P_INERTIA
       DIMENSION DOWN(40959,8), LABEL(40959,8), NP(40959,8)
       COMMON /TREE_INTG/ NP, RR15, LABEL, RR16, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, RR17, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
       DIMENSION D_X(3)
C*$* Padding ( RR17, RR16, RR15, RR14, RR13, RR12, RR11, RR10, RR9, RR8
C*$*& , RR7 )
       INTEGER*2 HI4, HI3, HI2, HI1
       INTEGER II3, II2, II1
       PARAMETER (HI4 = 5, HI3 = 3, HI2 = 2, II3 = 2, II2 = 1)
       PARAMETER (II1 = 3, HI1 = 1)
       REAL RR6, RR5, RR4, RR3, RR2, RR1
       REAL RR17 (86), RR16 (86), RR15 (86), RR14 (86), RR13 (86), RR12
     X    (86), RR11 (86), RR10 (86), RR9 (86), RR8 (86), RR7 (86)
!     square cell size (it's easy with a single bit value):
       S2 = AREST(NODE) * AREST(NODE)
!     calculate cell-to-particle squared hard-distance.
!     r2 <-- r * r
       R2 = 0.
       DO J=HI1,II1
        ! D_x is the position vector of the particle with respect to
        ! the observed cell:
        !
        RR5 = XP(I,J) - X_CELL(NODE,CELL,J)
        R2 = R2 + RR5 * RR5
        D_X(J) = RR5
       END DO
!     if this is not a good cell, then jump it out; otherwise, set this
!     one as a terminal cell:
       IF (NP(NODE,CELL) .GT. II2) THEN
        IF (THETA2 * R2 .LT. S2) THEN
         RETURN 
        END IF
        PERMIT(DOWN(NODE,CELL)) = .FALSE.
       ELSE
        IF (LABEL(NODE,CELL) .EQ. I) RETURN 
        !since we may not calculate particle's self gravity:
       END IF
!     count the good-cell entry:
       KWSC = KWSC + HI1
!     calculate softened cell-to-particle linear, cubic and squared
!     distances respectively:
       R2 = R2 + EPSILON2
       R = SQRT (R2)
       R3 = R2 * R
       RR1 = II2 / R2
       DO J=HI1,II1
        RR6 = D_X(J) * RR1
        D_X(J) = RR6
       END DO
       PHI0 = -CELLMASS(NODE,CELL) / R
       DO J=HI1,II1
        GNODE(KWSC,J) = PHI0 * D_X(J)
       END DO
       DIAG = 0.
       TRACE = 0.
       DO K=HI1,II1
        DIAG = DIAG + D_X(K) * P_INERTIA(NODE,CELL,K,K) * D_X(K)
        TRACE = TRACE + P_INERTIA(NODE,CELL,K,K)
       END DO
       TRIANG = 0.
       DO L=II2,II3
        RR2 = D_X(L)
        DO K=L+II2,II1
         TRIANG = TRIANG + D_X(K) * P_INERTIA(NODE,CELL,K,L) * RR2
        END DO
       END DO
       DOUB = DIAG + TRIANG * HI2
       PHI2 = -.5 * (DOUB * HI3 - TRACE / R2) / R
       RR3 = PHI2 * HI4
       RR4 = II2 / R3
       DO J=HI1,II1
        SINGL = 0.
        DO K=HI1,II1
         SINGL = SINGL + P_INERTIA(NODE,CELL,J,K) * D_X(K)
        END DO
        GNODE(KWSC,J) = GNODE(KWSC,J) + (SINGL * HI3 - D_X(J) * TRACE) *
     X     RR4 + RR3 * D_X(J)
       END DO
      END
 
 
# 407 "SPH_grav.F"
 
