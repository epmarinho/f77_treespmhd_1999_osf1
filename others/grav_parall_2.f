! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 13:16:04
!### KAP/Digital_UA_F detected intrinsic conflicts in this program.
!### Check the KAP/Digital_UA_F listing for details.
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:15:25
# 1 "SPH_grav.F"
 
! 
!                  TREE-GRAVITY SUBROUTINES
! 
 
 
      SUBROUTINE TREEPOTENTIAL
       INTEGER CELL, I
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
       COMMON /MOMENTA/ VMOMENT, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, U_OLD, U_DOT, U_DOT_VISC, LAMBDA, 
     X   VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, KSEED, LSEED, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, Q_ACC, P_ACC
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
       COMMON /TREE_INTG/ NP, LABEL, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
       INTEGER*2 HI2, HI1
       INTEGER II3, II2
       PARAMETER (HI2 = 0, II3 = 8, II2 = 0, HI1 = 1)
       EXTERNAL MPPTID, MPTEPA, MPTXPA, MPPFKD, MPOFRK, PKTREEPOTENTIAL_
       INTEGER MPPTID, MPPFKD, II9(0:10)
       CHARACTER AA4*30
       INTEGER II8, II1, NODE2, NODE1
       AUTOMATIC II8
       DATA II9(0)/0/ 
       DATA II9(1)/0/ 
       DATA II9(2)/30/ 
       DATA AA4/';SPH_grav.F;TREEPOTENTIAL;13;;'/ 
       EQUIVALENCE (AA4, II9(3))
       INTEGER II11, II10
       PARAMETER (II11 = 0, II10 = 500)
       II8 = MPPTID ()
       W_TOT = 0.
       W_S = 0.
       DO I=HI1,N
        IF (MAXNODE .GT. II10 .AND. MPPFKD () .EQ. II11) THEN
         CALL MPOFRK (PKTREEPOTENTIAL_,II9)
        ELSE
C!!!!! PARALLEL REGION IF (MAXNODE .GT. 500) SHARED (MAXNODE,PERMIT) 
C!!!!!& LOCAL (NODE1)
         CALL MPTEPA (II8)
C!!!!! PARALLEL DO 
         DO NODE1=HI1,MAXNODE
          PERMIT(NODE1) = .TRUE.
         END DO
C!!!!! END PARALLEL DO NOWAIT
         CALL MPTXPA (II8)
C!!!!! END PARALLEL REGION 
        END IF
        IF (MAXNODE .GT. II2) NODE = MAXNODE
        II1 = MAXNODE
        DO NODE2=HI1,II1
         NODE = NODE2
          !this descent-control breaks vectorization
         IF (PERMIT(NODE2)) THEN
          DO CELL=HI1,II3
           IF (NP(NODE,CELL) .GT. HI2) CALL POTENTIAL (I,CELL)
          END DO
         ELSE
          CALL FORBID
         END IF
        END DO
        NODE = MAX0 (II1, II2) + HI1
       END DO
       RETURN 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:16:04
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:15:25
      SUBROUTINE PKTREEPOTENTIAL_ ( MPPID, MPPNPR )
       INTEGER*2 HI11
       PARAMETER (HI11 = 1)
       INTEGER NODE3
       INTEGER MPPID
       INTEGER NNODE1
       PARAMETER (NNODE1 = 40959)
       AUTOMATIC II5
       INTEGER II5
       AUTOMATIC II7
       INTEGER II7
       AUTOMATIC NODE4
       INTEGER NODE4
       INTEGER MPPNPR
       INTEGER MAXNODE1
       AUTOMATIC II4
       INTEGER II4
       AUTOMATIC II6
       INTEGER II6
       LOGICAL PERMIT1(40959)
       COMMON /NODE/ NODE3, MAXNODE1
       COMMON /PERMIT/ PERMIT1
C!!!!! PARALLEL REGION IF (MAXNODE .GT. 500) SHARED (MAXNODE,PERMIT) 
C!!!!!& LOCAL (NODE1)
       INTEGER II12
       PARAMETER (II12 = 1)
       INTEGER*2 HI1
       INTEGER II1
       PARAMETER (HI1 = 1, II1 = 1)
       EXTERNAL MPPTID, MPTEPA, MPTXPA, MPPFKD, MPOFRK, 
     X   PKPKTREEPOTENTIAL__
       INTEGER MPPTID, MPPFKD, II13(0:11)
       CHARACTER AA1*36
       INTEGER II11
       AUTOMATIC II11
       DATA II13(0)/3/ 
       DATA II13(1)/0/ 
       DATA II13(2)/36/ 
       DATA AA1/';SPH_common.F;PKTREEPOTENTIAL_;514;;'/ 
       EQUIVALENCE (AA1, II13(3))
       INTEGER II15, II14
       PARAMETER (II15 = 0, II14 = 500)
       II11 = MPPTID ()
       II7 = (MAXNODE1 + MPPNPR - II1) / MPPNPR
       II5 = MIN (MAXNODE1, MPPID * II7 + II7)
C!!!!! PARALLEL DO 
       IF (II5 .GT. MPPID * II7 + II14 .AND. MPPFKD () .EQ. II15) THEN
        CALL MPOFRK (PKPKTREEPOTENTIAL__,II13,II7,MPPID,II5)
       ELSE
C!!!!! PARALLEL REGION IF (II5 .GT. MPPID * II7 + 500) SHARED (MPPID,II7
C!!!!!& ,II5,PERMIT1) LOCAL (NODE4)
        CALL MPTEPA (II11)
C!!!!! PARALLEL DO 
        DO NODE4=MPPID*II7+HI1,II5
         PERMIT1(NODE4) = .TRUE.
        END DO
C!!!!! END PARALLEL DO NOWAIT
        CALL MPTXPA (II11)
C!!!!! END PARALLEL REGION 
       END IF
C!!!!! END PARALLEL DO NOWAIT
C!!!!! END PARALLEL REGION 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:16:04
C!!!!! PARALLEL DO 
      SUBROUTINE PKPKTREEPOTENTIAL__ ( MPPID, MPPNPR, II71, MPPID1, II51
     X   )
       AUTOMATIC NODE41
       INTEGER NODE41
       INTEGER MPPID
       INTEGER II71
       AUTOMATIC II8
       INTEGER II8, MPPTID
       EXTERNAL MPPTID
       LOGICAL PERMIT11(40959)
       INTEGER*2 HI12
       PARAMETER (HI12 = 1)
       AUTOMATIC II2
       INTEGER II2
       INTEGER II51
       AUTOMATIC II9
       INTEGER II9
       AUTOMATIC II10
       INTEGER II10
       INTEGER*8 MPPFOD1(0:7,0:127)
       INTEGER MPPNPR
       INTEGER MPPID1
       AUTOMATIC II3
       INTEGER II3, MPPIOA
       EXTERNAL MPPIOA
       INTEGER MPPFOA1
       COMMON /MPPFOA/ MPPFOA1
       COMMON /MPPFOD/ MPPFOD1
       COMMON /PERMIT/ PERMIT11
       INTEGER II19, II18, II17, II16
       PARAMETER (II19 = 1, II18 = 2, II17 = 3, II16 = 0)
       II10 = MPPTID ()
       IF (MPPFOA1 .LT. II16) THEN
        MPPFOD1(MPPIOA (II17),II10) = %LOC (II51)
        MPPFOD1(MPPIOA (II18),II10) = %LOC (MPPID1)
        MPPFOD1(MPPIOA (II19),II10) = %LOC (II71)
       END IF
C!!!!! PARALLEL REGION IF (II5 .GT. MPPID * II7 + 500) SHARED (MPPID,II7
C!!!!!& ,II5,PERMIT1) LOCAL (NODE4)
       II8 = II51 - (MPPID1 * II71 + HI12) + II19
       II9 = (II8 + MPPNPR - II19) / MPPNPR
       II2 = MPPID1 * II71 + HI12 + MPPID * II9
       II3 = MIN (II51, II2 + (II9 - II19))
C!!!!! PARALLEL DO 
       DO NODE41=II2,II3
        PERMIT11(NODE41) = .TRUE.
       END DO
C!!!!! END PARALLEL DO NOWAIT
C!!!!! END PARALLEL REGION 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:16:04
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:15:25
 
 
      SUBROUTINE POTENTIAL ( I, CELL )
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
       COMMON /MOMENTA/ VMOMENT, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, U_OLD, U_DOT, U_DOT_VISC, LAMBDA, 
     X   VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, KSEED, LSEED, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, Q_ACC, P_ACC
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
       COMMON /TREE_INTG/ NP, LABEL, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
       DIMENSION AUX(40960,3)
       COMMON /AUXILIAR/ AUX
       EQUIVALENCE (AUX, PS)
       DIMENSION PS(40960,3)
       INTEGER*2 HI3, HI2, HI1
       INTEGER II3, II2, II1
       PARAMETER (HI3 = 3, HI2 = 2, II3 = 2, II2 = 1, II1 = 3)
       PARAMETER (HI1 = 1)
       REAL RR3, RR2, RR1
       INTEGER*2 HI6, HI5, HI4
       INTEGER II6, II5, II4
       PARAMETER (HI6 = 3, HI5 = 2, II6 = 2, II5 = 1, II4 = 3)
       PARAMETER (HI4 = 1)
! 
!      Part I. Searching for a well separated cell (good cell)
! 
!      square up the box size:
       S2 = AREST(NODE)
       S2 = S2 * S2
!      calculate cell-to-particle actual squared distance.
!      r2 <-- r * r
       R2 = 0.
       DO J=HI4,II4
        RR2 = XP(I,J) - X_CELL(NODE,CELL,J)
        R2 = R2 + RR2 * RR2
        PS(HI4,J) = RR2
       END DO
!      if it is not a good cell then go out
       IF (NP(NODE,CELL) .GT. II5) THEN
        IF (THETA2 * R2 .LT. S2) THEN
         RETURN 
        END IF
        PERMIT(DOWN(NODE,CELL)) = .FALSE.
       ELSE
        IF (LABEL(NODE,CELL) .EQ. I) RETURN 
       END IF
! 
!      Part II. Perform cell contribution to collective potential 
! 
 
 
 
 
 
 
 
 
 
 
!      square softened distance:
       R2 = R2 + EPSILON2
!      softened distance:
       R = SQRT (R2)
!      cubic softened distance:
       R3 = R2 * R
       RR1 = II5 / R2
!      calculate pole-vector:
       DO J=HI4,II4
        RR3 = PS(HI4,J) * RR1
        PS(HI4,J) = RR3
       END DO
!      --------------------------------
!      perform monopole approximation
!      ------------------------------
       PHI0 = -CELLMASS(NODE,CELL) / R
!      --------------------------------
!      perform quadrupole approximation
!      --------------------------------
!      diagonal sumation:
       DIAG = 0.
       TRACE = 0.
       DO K=HI4,II4
        DIAG = DIAG + PS(HI4,K) * P_INERTIA(NODE,CELL,K,K) * PS(HI4,K)
        TRACE = TRACE + P_INERTIA(NODE,CELL,K,K)
       END DO
!      triangular sumation: quadrupole is symmetric!
       TRIANG = 0.
       DO L=II5,II6
        DO K=L+II5,II4
         TRIANG = TRIANG + PS(HI4,K) * P_INERTIA(NODE,CELL,K,L) * PS(HI4
     X     ,L)
        END DO
       END DO
!      doub gets the double internal-product:    ^p . ^Q^ . ^p
       DOUB = DIAG + TRIANG * HI5
!      perform quadrupole contribution:
       PHI2 = -.5 * (DOUB * HI6 - TRACE / R2) / R
!      calculate potential contribution
       PHI = PHI0 + PHI2
       E_PAIR = MASS(I) * PHI
!      add "E_pair" to total potential-energy "W_tot"
       W_TOT = W_TOT + .5 * E_PAIR
!      calculate soft-virial contribution
       E_S = .5 * MASS(I) * CELLMASS(NODE,CELL) * EPSILON2 / R3
!      add "E_s" to total soft-virial "W_s"
       W_S = W_S + E_S
       RETURN 
      END
! KAP/Digital_UA_F 4.0 k3011126 980529  17-Jul-2000 13:16:04
!### KAP/Digital_UA_F detected intrinsic conflicts in this program.
!### Check the KAP/Digital_UA_F listing for details.
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:15:25
 
 
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
       COMMON /MOMENTA/ VMOMENT, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, U_OLD, U_DOT, U_DOT_VISC, LAMBDA, 
     X   VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, KSEED, LSEED, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, Q_ACC, P_ACC
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
       COMMON /TREE_INTG/ NP, LABEL, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
       INTEGER II17, II16, II15, II14, II13, II12, II11
       INTEGER*2 HI2, HI1
       PARAMETER (II17 = 2, II16 = 4, II15 = 3, II14 = 2730, II13 = 1)
       PARAMETER (HI2 = 0, II12 = 8, II11 = 0, HI1 = 1)
       EXTERNAL MPPTID, MPTEPA, MPTXPA, MPPFKD, MPOFRK, PKTREEGRAVITY_
       INTEGER MPPTID, MPPFKD, II23(0:10)
       CHARACTER AA8*29
       INTEGER II22, II10, II9, II8, II7, II6, II5, II4, II3, II2, II1, 
     X   NODE2, NODE1
       REAL RR2, RR1
       LOGICAL LL1
       AUTOMATIC II22
       DATA II23(0)/0/ 
       DATA II23(1)/0/ 
       DATA II23(2)/29/ 
       DATA AA8/';SPH_grav.F;TREEGRAVITY;129;;'/ 
       EQUIVALENCE (AA8, II23(3))
       INTEGER II25, II24
       PARAMETER (II25 = 0, II24 = 500)
       II22 = MPPTID ()
      !dimension xa(dim),ga(dim)
       PRINT *, 'treegravity'
       DO IL=HI1,N_P(TIMEBIN)
        I = P_LIST(IL,TIMEBIN)
        IF (MAXNODE .GT. II24 .AND. MPPFKD () .EQ. II25) THEN
         CALL MPOFRK (PKTREEGRAVITY_,II23)
        ELSE
C!!!!! PARALLEL REGION IF (MAXNODE .GT. 500) SHARED (MAXNODE,PERMIT) 
C!!!!!& LOCAL (NODE1)
         CALL MPTEPA (II22)
C!!!!! PARALLEL DO 
         DO NODE1=HI1,MAXNODE
          PERMIT(NODE1) = .TRUE.
         END DO
C!!!!! END PARALLEL DO NOWAIT
         CALL MPTXPA (II22)
C!!!!! END PARALLEL REGION 
        END IF
        IF (MAXNODE .GT. II11) NODE = MAXNODE
        KWSC = II11
        ! Performing tree-walk
        II1 = MAXNODE
        DO NODE2=HI1,II1
         NODE = NODE2
          ! vectorization is broken here:
         IF (PERMIT(NODE2)) THEN
          DO CELL=HI1,II12
           IF (NP(NODE,CELL) .GT. HI2) CALL GRAVITY (I,CELL)
          END DO
         ELSE
          CALL FORBID
         END IF
        END DO
        NODE = MAX0 (II1, II11) + HI1
        II5 = MOD (KWSC - II13, II14) + II13
        !-----------------------------
        !do j=1,dim
        !  xa(j)=xp(i,j)
        !end do
        !call phi_galax(xa,ga,epsilon)
        ! Now, the particle has a list
        ! of kwsc cell'contributions:
        !-----------------------------
        DO J=HI1,II15
         G(I,J) = 0.
        END DO
        II3 = II13
        II2 = II5
        DO II6=II13,KWSC,II14
         II4 = II3 + II2 - II13
         II7 = (II4 - II3 + II13) / II16
         LL1 = II7 .GT. II11
         II10 = II7 * II16 + II3
        !-----------------------------
        !do j=1,dim
        !  xa(j)=xp(i,j)
        !end do
        !call phi_galax(xa,ga,epsilon)
        ! Now, the particle has a list
        ! of kwsc cell'contributions:
        !-----------------------------
         DO J=HI1,II15
          IF (LL1) THEN
           RR1 = G(I,J)
           DO L=II3,II4-II15,II16
            RR1 = RR1 + GNODE(L,J)
            RR1 = RR1 + GNODE(L+II13,J)
            RR1 = RR1 + GNODE(L+II17,J)
            RR1 = RR1 + GNODE(L+II15,J)
            !g(i,j)=g(i,j)+ga(j)
           END DO
           G(I,J) = RR1
          END IF
          L = II10
          II8 = L
          II9 = II4 - II8 + II13
          IF (II9 .GT. II11) THEN
           RR2 = G(I,J)
           DO L=II8,II4,II13
            RR2 = RR2 + GNODE(L,J)
           END DO
           G(I,J) = RR2
          END IF
         END DO
         II3 = II4 + II13
         II2 = II14
        END DO
        !-----------------------------
       END DO
       PRINT *, 'done!'
       RETURN 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:16:04
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:15:25
      SUBROUTINE PKTREEGRAVITY_ ( MPPID, MPPNPR )
       AUTOMATIC II20
       INTEGER II20
       INTEGER MPPID
       AUTOMATIC NODE4
       INTEGER NODE4
       INTEGER NODE3
       AUTOMATIC II21
       INTEGER II21
       AUTOMATIC II18
       INTEGER II18
       INTEGER*2 HI11
       PARAMETER (HI11 = 1)
       INTEGER MAXNODE1
       AUTOMATIC II19
       INTEGER II19
       INTEGER MPPNPR
       INTEGER NNODE1
       PARAMETER (NNODE1 = 40959)
       LOGICAL PERMIT1(40959)
       COMMON /NODE/ NODE3, MAXNODE1
       COMMON /PERMIT/ PERMIT1
C!!!!! PARALLEL REGION IF (MAXNODE .GT. 500) SHARED (MAXNODE,PERMIT) 
C!!!!!& LOCAL (NODE1)
       INTEGER II26
       PARAMETER (II26 = 1)
       INTEGER*2 HI1
       INTEGER II1
       PARAMETER (HI1 = 1, II1 = 1)
       EXTERNAL MPPTID, MPTEPA, MPTXPA, MPPFKD, MPOFRK, 
     X   PKPKTREEGRAVITY__
       INTEGER MPPTID, MPPFKD, II8(0:11)
       CHARACTER AA1*34
       INTEGER II7
       AUTOMATIC II7
       DATA II8(0)/3/ 
       DATA II8(1)/0/ 
       DATA II8(2)/34/ 
       DATA AA1/';SPH_common.F;PKTREEGRAVITY_;576;;'/ 
       EQUIVALENCE (AA1, II8(3))
       INTEGER II10, II9
       PARAMETER (II10 = 0, II9 = 500)
       II7 = MPPTID ()
       II21 = (MAXNODE1 + MPPNPR - II1) / MPPNPR
       II19 = MIN (MAXNODE1, MPPID * II21 + II21)
C!!!!! PARALLEL DO 
       IF (II19 .GT. MPPID * II21 + II9 .AND. MPPFKD () .EQ. II10) THEN
        CALL MPOFRK (PKPKTREEGRAVITY__,II8,MPPID,II19,II21)
       ELSE
C!!!!! PARALLEL REGION IF (II19 .GT. MPPID * II21 + 500) SHARED (MPPID,
C!!!!!& II21,II19,PERMIT1) LOCAL (NODE4)
        CALL MPTEPA (II7)
C!!!!! PARALLEL DO 
        DO NODE4=MPPID*II21+HI1,II19
         PERMIT1(NODE4) = .TRUE.
        END DO
C!!!!! END PARALLEL DO NOWAIT
        CALL MPTXPA (II7)
C!!!!! END PARALLEL REGION 
       END IF
C!!!!! END PARALLEL DO NOWAIT
C!!!!! END PARALLEL REGION 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:16:04
C!!!!! PARALLEL DO 
      SUBROUTINE PKPKTREEGRAVITY__ ( MPPID, MPPNPR, MPPID1, II191, II211
     X   )
       AUTOMATIC NODE41
       INTEGER NODE41
       INTEGER MPPID
       INTEGER II211
       AUTOMATIC II4
       INTEGER II4, MPPTID
       EXTERNAL MPPTID
       LOGICAL PERMIT11(40959)
       INTEGER*2 HI12
       PARAMETER (HI12 = 1)
       AUTOMATIC II2
       INTEGER II2
       INTEGER II191
       AUTOMATIC II5
       INTEGER II5
       AUTOMATIC II6
       INTEGER II6
       INTEGER*8 MPPFOD1(0:7,0:127)
       INTEGER MPPNPR
       INTEGER MPPID1
       AUTOMATIC II3
       INTEGER II3, MPPIOA
       EXTERNAL MPPIOA
       INTEGER MPPFOA1
       COMMON /MPPFOA/ MPPFOA1
       COMMON /MPPFOD/ MPPFOD1
       COMMON /PERMIT/ PERMIT11
       INTEGER II14, II13, II12, II11
       PARAMETER (II14 = 1, II13 = 2, II12 = 3, II11 = 0)
       II6 = MPPTID ()
       IF (MPPFOA1 .LT. II11) THEN
        MPPFOD1(MPPIOA (II12),II6) = %LOC (II211)
        MPPFOD1(MPPIOA (II13),II6) = %LOC (II191)
        MPPFOD1(MPPIOA (II14),II6) = %LOC (MPPID1)
       END IF
C!!!!! PARALLEL REGION IF (II19 .GT. MPPID * II21 + 500) SHARED (MPPID,
C!!!!!& II21,II19,PERMIT1) LOCAL (NODE4)
       II4 = II191 - (MPPID1 * II211 + HI12) + II14
       II5 = (II4 + MPPNPR - II14) / MPPNPR
       II2 = MPPID1 * II211 + HI12 + MPPID * II5
       II3 = MIN (II191, II2 + (II5 - II14))
C!!!!! PARALLEL DO 
       DO NODE41=II2,II3
        PERMIT11(NODE41) = .TRUE.
       END DO
C!!!!! END PARALLEL DO NOWAIT
C!!!!! END PARALLEL REGION 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:16:04
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 13:15:25
 
 
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
       COMMON /MOMENTA/ VMOMENT, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, U_OLD, U_DOT, U_DOT_VISC, LAMBDA, 
     X   VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, KSEED, LSEED, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, Q_ACC, P_ACC
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
       COMMON /TREE_INTG/ NP, LABEL, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
       DIMENSION AUX(40960,3)
       COMMON /AUXILIAR/ AUX
       DIMENSION D_X(40960,3)
       EQUIVALENCE (AUX, D_X)
       INTEGER*2 HI4, HI3, HI2, HI1
       INTEGER II3, II2, II1
       PARAMETER (HI4 = 5, HI3 = 3, HI2 = 2, II3 = 2, II2 = 1)
       PARAMETER (II1 = 3, HI1 = 1)
       REAL RR5, RR4, RR3, RR2, RR1
       INTEGER*2 HI8, HI7, HI6, HI5
       INTEGER II6, II5, II4
       PARAMETER (HI8 = 5, HI7 = 3, HI6 = 2, II6 = 2, II5 = 1)
       PARAMETER (II4 = 3, HI5 = 1)
**************************************************************
*     performs a tree-descent by searching for well-separated
*     cells in the tree
**************************************************************
!     square cell size (it's easy with a single bit value):
       S2 = AREST(NODE) * AREST(NODE)
!     calculate cell-to-particle squared hard-distance.
!     r2 <-- r * r
       R2 = 0.
       DO J=HI5,II4
        ! D_x is the position vector of the particle with respect to
        ! the observed cell:
        !
        RR4 = XP(I,J) - X_CELL(NODE,CELL,J)
 
 
 
 
 
 
 
 
 
 
 
 
 
        R2 = R2 + RR4 * RR4
        D_X(HI5,J) = RR4
       END DO
!     if this is not a good cell, then jump it out; otherwise, set this
!     one as a terminal cell:
       IF (NP(NODE,CELL) .GT. II5) THEN
        !** itheta=-2*int(s(i)/sqrt(r2))
        !** theta_tmp=theta2*(2**itheta)
 
        !** print*,'theta_tmp=',theta_tmp,' itheta=',itheta
 
        !** if((theta_tmp*r2).lt.s2)then
        IF (THETA2 * R2 .LT. S2) THEN
         RETURN 
        END IF
        PERMIT(DOWN(NODE,CELL)) = .FALSE.
       ELSE
        IF (LABEL(NODE,CELL) .EQ. I) RETURN 
        !since we may not calculate particle's self gravity:
       END IF
**************************************************************
!     count the good-cell entry:
       KWSC = KWSC + HI5
!     calculate softened cell-to-particle linear, cubic and squared
!     distances respectively:
 
 
 
 
 
 
 
 
 
 
 
 
       R2 = R2 + EPSILON2
       R = SQRT (R2)
       RR1 = II5 / R2
!      calculate "pole-vector" (softened unit vector):
       DO J=HI5,II4
        RR5 = D_X(HI5,J) * RR1
        D_X(HI5,J) = RR5
       END DO
 
! **************************
! perform cell contributions
! **************************
 
!     ******************************
!     perform monopole approximation
!     ******************************
 
!     potential:
       PHI0 = -CELLMASS(NODE,CELL) / R
      !for testing:
      !*** if(phi0.eq.0)then
            !*** print*,'gravity: error: cellmass is zero'
            !*** stop
      !*** endif
 
!      monopole-gravity-field:
       DO J=HI5,II4
        GNODE(KWSC,J) = PHI0 * D_X(HI5,J)
       END DO
 
!     ********************************
!     perform quadrupole approximation
!     ********************************
 
!     diagonal sumation:
       DIAG = 0.
       TRACE = 0.
       DO K=HI5,II4
        DIAG = DIAG + D_X(HI5,K) * P_INERTIA(NODE,CELL,K,K) * D_X(HI5,K)
        TRACE = TRACE + P_INERTIA(NODE,CELL,K,K)
       END DO
!      triangular sumation: quadrupole is symetric!
       TRIANG = 0.
       DO L=II5,II6
        DO K=L+II5,II4
         TRIANG = TRIANG + D_X(HI5,K) * P_INERTIA(NODE,CELL,K,L) * D_X(
     X     HI5,L)
        END DO
       END DO
!      doub gets double the internal-product:    ^p . ^q^ . ^p
       DOUB = DIAG + TRIANG * HI6
!      perform quadrupole cell contribution to potential:
       PHI2 = -.5 * (DOUB * HI7 - TRACE / R2) / R
!      calculate potential contribution
!      total = monopole + quadrupole
       RR2 = PHI2 * HI8
       RR3 = II5 / R
 
       DO J=HI5,II4
!        singl gets the single internal-product:    ^q^ . ^p
        SINGL = 0.
        DO K=HI5,II4
         SINGL = SINGL + P_INERTIA(NODE,CELL,J,K) * D_X(HI5,K)
        END DO
        GNODE(KWSC,J) = GNODE(KWSC,J) + (SINGL * HI7 - D_X(HI5,J) * 
     X    TRACE) * RR3 + RR2 * D_X(HI5,J)
       END DO
       RETURN 
      END
 
 
# 407 "SPH_grav.F"
 
