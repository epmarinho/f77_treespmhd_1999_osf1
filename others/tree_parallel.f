C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
# 1 "SPH_tree.F"
! 
!
!      TREE CONSTRUCTION AND TREE-ACCESS SUBROUTINES:
!
! 
!
      SUBROUTINE SETBOX
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
       COMMON /MOMENTA/ VMOMENT, RR1, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, RR2, U_OLD, U_DOT, RR3, U_DOT_VISC
     X   , LAMBDA, RR4, VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, RR5, KSEED, LSEED, RR6, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, RR7, Q_ACC, RR8, P_ACC
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
       COMMON /TREE_INTG/ NP, RR9, LABEL, RR10, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, RR11, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
       REAL AUX, VERTIX
 
 
 
       REAL SCALE, S1BIT, AR, C, A
       DIMENSION AUX(40960,3)
 
       COMMON /AUX/ AUX
       EQUIVALENCE (AUX, VERTIX)
       DIMENSION VERTIX(1,3)
C*$* Padding ( RR11, RR10, RR9, RR8, RR7, RR6, RR5, RR4, RR3, RR2, RR1 )
       INTEGER*2 HI3, HI2, HI1
       INTEGER II1
       PARAMETER (HI3 = 0, HI2 = 2, II1 = 3, HI1 = 1)
       EXTERNAL MPPTID, MPTEPA, MPTXPA, MPPFKD, MPOFRK, PKSETBOX_
       INTEGER MPPTID, MPPFKD, II9(0:8)
       CHARACTER AA4*23
       INTEGER II8
       REAL RR11 (86), RR10 (86), RR9 (86), RR8 (86), RR7 (86), RR6 (86)
     X   , RR5 (86), RR4 (86), RR3 (86), RR2 (86), RR1 (86)
       AUTOMATIC II8
       DATA II9(0)/1/ 
       DATA II9(1)/0/ 
       DATA II9(2)/23/ 
       DATA AA4/';SPH_tree.F;SETBOX;69;;'/ 
       EQUIVALENCE (AA4, II9(3))
       INTEGER II12, II10
       PARAMETER (II12 = 0, II10 = 166)
       II8 = MPPTID ()
!
!
!      * toma o vertice inferior do triedro circunscrito no sistema:
!          vertix(j)=min(xp(i,j))
!      * determina a maior aresta do cubo inscrito no triedro:
!          arest = max(xp(i,j))-min(xp(i,j))
!
       DO J=HI1,II1
        VERTIX(HI1,J) = XP(HI1,J)
       END DO
       C = XP(HI1,HI1)
       A = C
       DO J=HI1,II1
        DO I=HI2,N
         IF (XP(I,J) .GT. C) THEN
          C = XP(I,J)
         ELSE
          VERTIX(HI1,J) = MIN (VERTIX(HI1,J), XP(I,J))
         END IF
        END DO
       END DO
       DO J=HI1,II1
        A = MIN (A, VERTIX(HI1,J))
       END DO
       AR = C - A
!
!      * determina o cubo circunscrito cujas arestas sao potencias
!        inteiras de 2:
!
       S1BIT = 1.
       IF (AR .GT. 1.) THEN
        DO WHILE ( S1BIT .LT. AR )
         S1BIT = S1BIT * 2.
        END DO
       ELSE
        AR = AR * HI2
        DO WHILE ( S1BIT .GT. AR )
         S1BIT = 0.5 * S1BIT
        END DO
       END IF
       IF (S1BIT .NE. HI3) THEN
        SCALE = 1D0 / S1BIT
        AREST(HI1) = 0.5 * S1BIT
       ELSE
        PRINT *, '%err: particles were distributed into a singularity!!'
        STOP 
       END IF
       IF (N .GT. II10 .AND. MPPFKD () .EQ. II12) THEN
        CALL MPOFRK (PKSETBOX_,II9,SCALE)
       ELSE
C!!!!! PARALLEL REGION IF (N .GT. 166) SHARED (N,BOX,SCALE,XP,VERTIX) 
C!!!!!& LOCAL (I,J)
        CALL MPTEPA (II8)
!
!      renormaliza as coordenadas das particulas para um cubo
!      unitario:
!
        DO J=HI1,II1
C!!!!! PARALLEL DO 
         DO I=HI1,N
 
 
 
 
          BOX(I,J) = SCALE * (XP(I,J) - VERTIX(HI1,J))
 
         END DO
C!!!!! END PARALLEL DO NOWAIT
C!!!!! BARRIER
        END DO
        CALL MPTXPA (II8)
C!!!!! END PARALLEL REGION 
       END IF
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
      SUBROUTINE PKSETBOX_ ( MPPID, MPPNPR, SCALE1 )
       EQUIVALENCE (VERTIX1(1,1), RR12(1))
       EQUIVALENCE (AUX1(1,1), RR12(1))
       AUTOMATIC II3
       INTEGER II3, MPPTID
       EXTERNAL MPPTID
       INTEGER MPPID
       INTEGER NN1
       PARAMETER (NN1 = 40960)
       REAL RR12(122880)
       INTEGER*8 MPPFOD1(0:7,0:127)
       INTEGER II11
       PARAMETER (II11 = 3)
       INTEGER MPPFOA1
       REAL SCALE1
       AUTOMATIC II4
       INTEGER II4
       AUTOMATIC I1
       INTEGER I1
       AUTOMATIC II7
       INTEGER II7
       AUTOMATIC II2
       INTEGER II2, MPPIOA
       EXTERNAL MPPIOA
       CHARACTER AA3*23
       INTEGER DIM1
       PARAMETER (DIM1 = 3)
       INTEGER MPPNPR
       INTEGER N1
       INTEGER II6(0:8)
       INTEGER*2 HI11
       PARAMETER (HI11 = 1)
       AUTOMATIC II5
       INTEGER II5
       AUTOMATIC J1
       INTEGER J1
       EXTERNAL MPSBAR
       COMMON /MPPFOD/ MPPFOD1
       COMMON /TOTAL_NUMBER/ N1
       REAL XP1(40960,3)
       REAL BOX1(40960,3)
       REAL VERTIX1(1,3)
       REAL AUX1(40960,3)
       COMMON /MPPFOA/ MPPFOA1
       COMMON /BOXY_COORDINATES/ BOX1
       COMMON /AUX/ AUX1
       COMMON /POSITIONS/ XP1
       DATA II6(0)/0/ 
       DATA II6(1)/0/ 
       DATA II6(2)/23/ 
       DATA AA3/';SPH_tree.F;SETBOX;76;;'/ 
       EQUIVALENCE (AA3, II6(3))
       INTEGER II14, II13
       PARAMETER (II14 = 1, II13 = 0)
       II7 = MPPTID ()
       IF (MPPFOA1 .LT. II13) MPPFOD1(MPPIOA (II14),II7) = %LOC (SCALE1)
C!!!!! PARALLEL REGION IF (N .GT. 166) SHARED (N,BOX,SCALE,XP,VERTIX) 
C!!!!!& LOCAL (I,J)
!
!      renormaliza as coordenadas das particulas para um cubo
!      unitario:
!
       DO J1=HI11,II11
        II4 = N1 - HI11 + II14
        II5 = (II4 + MPPNPR - II14) / MPPNPR
        II2 = HI11 + MPPID * II5
        II3 = MIN (N1, II2 + (II5 - II14))
C!!!!! PARALLEL DO 
        DO I1=II2,II3
 
 
 
 
         BOX1(I1,J1) = SCALE1 * (XP1(I1,J1) - VERTIX1(HI11,J1))
 
        END DO
C!!!!! END PARALLEL DO NOWAIT
C!!!!! BARRIER
        CALL MPSBAR (II7,II6)
       END DO
C!!!!! END PARALLEL REGION 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
!
!
      SUBROUTINE MAKETREE
!      this procedure generates a tree data-structure for a given 
!      set of n mass-points in 3-dim space. the method consists on 
!      updating the particles in all corresponding nodes in each level
!      at once, which implies in a one descenting construction, 
!      discarding recursive procedures. the philosophy of tree-
!      descendent construction has the principal advantage on its easy
!      adaptation for vector machines.
!
!    data configuration:
!      nn=n
!            nnode = n - 1
!            maxoctants = eight * nnode
!
!
!
!
       INTEGER I, I1, NSAVE
       REAL VOLUME, DENS_MEAN
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
       COMMON /MOMENTA/ VMOMENT, RR1, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, RR2, U_OLD, U_DOT, RR3, U_DOT_VISC
     X   , LAMBDA, RR4, VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, RR5, KSEED, LSEED, RR6, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, RR7, Q_ACC, RR8, P_ACC
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
       COMMON /TREE_INTG/ NP, RR9, LABEL, RR10, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, RR11, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
       COMMON /INDICES/ I, I1
C*$* Padding ( RR11, RR10, RR9, RR8, RR7, RR6, RR5, RR4, RR3, RR2, RR1 )
       INTEGER II4, II3, II2, II1
       INTEGER*2 HI2, HI1
       PARAMETER (II4 = 2, II3 = 1, HI2 = 2, II2 = 4, II1 = 0)
       PARAMETER (HI1 = 1)
       EXTERNAL MPPTID, MPTEPA, MPTXPA, MPPFKD, MPOFRK, PKMAKETREE_
       INTEGER MPPTID, MPPFKD, II16(0:9)
       CHARACTER AA6*26
       INTEGER II15, I2
       REAL RR11 (86), RR10 (86), RR9 (86), RR8 (86), RR7 (86), RR6 (86)
     X   , RR5 (86), RR4 (86), RR3 (86), RR2 (86), RR1 (86)
       AUTOMATIC II15
       DATA II16(0)/0/ 
       DATA II16(1)/0/ 
       DATA II16(2)/26/ 
       DATA AA6/';SPH_tree.F;MAKETREE;112;;'/ 
       EQUIVALENCE (AA6, II16(3))
       INTEGER II17(0:9)
       CHARACTER AA7*26
       EXTERNAL PKMAKETREE_1
       DATA II17(0)/1/ 
       DATA II17(1)/0/ 
       DATA II17(2)/26/ 
       DATA AA7/';SPH_tree.F;MAKETREE;177;;'/ 
       EQUIVALENCE (AA7, II17(3))
       INTEGER II20, II19, II18
       PARAMETER (II20 = 250, II19 = 0, II18 = 166)
       II15 = MPPTID ()
!
!
!      inicializacao
!
 
 
 
 
       CALL SETBOX
       NSAVE = N
       IF (N .GT. II18 .AND. MPPFKD () .EQ. II19) THEN
        CALL MPOFRK (PKMAKETREE_,II16)
       ELSE
C!!!!! PARALLEL REGION IF (N .GT. 166) SHARED (N,T) LOCAL (I2)
        CALL MPTEPA (II15)
C!!!!! PARALLEL DO 
        DO I2=HI1,N
         T(I2,HI1) = II1
         T(I2,II2) = II1
         T(I2,HI2) = I2
        END DO
C!!!!! END PARALLEL DO NOWAIT
        CALL MPTXPA (II15)
C!!!!! END PARALLEL REGION 
       END IF
       I = MAX0 (N, II1) + HI1
 
 
 
 
 
 
 
 
 
       NODE = II3
       NEXT = II4
!
!      while any particle is non-leaf yet, continue building tree:
       DO WHILE ( N .GT. II3 )
        CALL FINDOCTANTS (II3,N)
        I = II3
        DO WHILE ( I .LE. N )
            !make node by sharing cubic-cells
            !with respective particles:
         CALL MAKENODE
 
 
 
 
 
 
!            link node to the tree:
         CALL LINKTREE
        END DO
        CALL NOLEAVES
       END DO
!
!      finalization of tree-construction:
       MAXNODE = NODE - HI1
       N = NSAVE
 
 
 
!
 
       IF (FIRST) THEN
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
        IF (NSAVE .GT. II20 .AND. MPPFKD () .EQ. II19) THEN
         CALL MPOFRK (PKMAKETREE_1,II17,NSAVE)
        ELSE
C!!!!! PARALLEL REGION IF (NSAVE .GT. 250) SHARED (NSAVE,H,S) LOCAL (I2)
         CALL MPTEPA (II15)
C!!!!! PARALLEL DO 
         DO I2=HI1,NSAVE
          H(I2) = 0.3 * S(I2)
         END DO
C!!!!! END PARALLEL DO NOWAIT
         CALL MPTXPA (II15)
C!!!!! END PARALLEL REGION 
        END IF
        I = MAX0 (NSAVE, II1) + HI1
 
        FIRST = .FALSE.
       END IF
 
 
 
 
 
 
 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
      SUBROUTINE PKMAKETREE_1 ( MPPID, MPPNPR, NSAVE1 )
       INTEGER*2 HI12
       PARAMETER (HI12 = 1)
       INTEGER NN2
       PARAMETER (NN2 = 40960)
       AUTOMATIC II9
       INTEGER II9, MPPTID
       EXTERNAL MPPTID
       AUTOMATIC II13
       INTEGER II13
       INTEGER*8 MPPFOD1(0:7,0:127)
       AUTOMATIC I4
       INTEGER I4
       INTEGER NSAVE1
       INTEGER MPPID, MPPIOA
       EXTERNAL MPPIOA
       AUTOMATIC II12
       INTEGER II12
       INTEGER MPPNPR
       AUTOMATIC II10
       INTEGER II10
       AUTOMATIC II14
       INTEGER II14
       INTEGER MPPFOA1
       REAL S1(40960)
       REAL H1(40960)
       COMMON /MPPFOA/ MPPFOA1
       COMMON /MPPFOD/ MPPFOD1
       REAL H_DOT1(40960)
       COMMON /PARTICLE_SIZE/ S1
       COMMON /SMOOTHING_LENGTHS/ H1, H_DOT1
       INTEGER II23, II22
       PARAMETER (II23 = 1, II22 = 0)
       II14 = MPPTID ()
       IF (MPPFOA1 .LT. II22) MPPFOD1(MPPIOA (II23),II14) = %LOC (NSAVE1
     X   )
C!!!!! PARALLEL REGION IF (NSAVE .GT. 250) SHARED (NSAVE,H,S) LOCAL (I2)
       II12 = NSAVE1 - HI12 + II23
       II13 = (II12 + MPPNPR - II23) / MPPNPR
       II9 = HI12 + MPPID * II13
       II10 = MIN (NSAVE1, II9 + (II13 - II23))
C!!!!! PARALLEL DO 
       DO I4=II9,II10
        H1(I4) = 0.3 * S1(I4)
       END DO
C!!!!! END PARALLEL DO NOWAIT
C!!!!! END PARALLEL REGION 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
      SUBROUTINE PKMAKETREE_ ( MPPID, MPPNPR )
       INTEGER*2 HI21
       PARAMETER (HI21 = 2)
       AUTOMATIC II8
       INTEGER II8
       INTEGER II21
       PARAMETER (II21 = 4)
       INTEGER FOUR1
       PARAMETER (FOUR1 = 4)
       AUTOMATIC II7
       INTEGER II7
       INTEGER II11
       PARAMETER (II11 = 0)
       INTEGER NN1
       PARAMETER (NN1 = 40960)
       AUTOMATIC II6
       INTEGER II6
       INTEGER*2 HI11
       PARAMETER (HI11 = 1)
       INTEGER N1
       AUTOMATIC II5
       INTEGER II5
       AUTOMATIC I3
       INTEGER I3
       INTEGER MPPID
       INTEGER MPPNPR
       COMMON /TOTAL_NUMBER/ N1
       INTEGER T3(40960,4)
       COMMON /T/ T3
C!!!!! PARALLEL REGION IF (N .GT. 166) SHARED (N,T) LOCAL (I2)
       INTEGER II24
       PARAMETER (II24 = 1)
       II7 = N1 - HI11 + II24
       II8 = (II7 + MPPNPR - II24) / MPPNPR
       II5 = HI11 + MPPID * II8
       II6 = MIN (N1, II5 + (II8 - II24))
C!!!!! PARALLEL DO 
       DO I3=II5,II6
        T3(I3,HI11) = II11
        T3(I3,II21) = II11
        T3(I3,HI21) = I3
       END DO
C!!!!! END PARALLEL DO NOWAIT
C!!!!! END PARALLEL REGION 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
 
 
 
      SUBROUTINE MAKENODE
       LOGICAL CHANGING
       INTEGER CELL, IP, I1, I, J
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
 
       COMMON /INDICES/ I, I1
       COMMON /NON_SINGULARITY/ ISINGULARITY
C*$* Padding ( RR13, RR12, RR11, RR10, RR9, RR8, RR7, RR6, RR5, RR4, RR3
C*$*&  )
       INTEGER II3, II2, II1
       INTEGER*2 HI3, HI2, HI1
       PARAMETER (II3 = 1, II2 = 3, HI3 = 3, HI2 = 1, HI1 = 2)
       PARAMETER (II1 = 0)
       SAVE CHANGING
       REAL RR2, RR1
       REAL RR13 (86), RR12 (86), RR11 (86), RR10 (86), RR9 (86), RR8 (
     X   86), RR7 (86), RR6 (86), RR5 (86), RR4 (86), RR3 (86)
       ISINGULARITY = II1
       CALL CLR
       I1 = I
       DO WHILE ( .NOT.NEXTNODE )
        IP = T(I,HI1)
        CELL = T(I,HI3) + HI2
        COUNT(CELL) = COUNT(CELL) + HI2
        CELLMASS(NODE,CELL) = CELLMASS(NODE,CELL) + MASS(IP)
        LST(COUNT(CELL),CELL) = IP
        RR1 = MASS(IP)
!         calculate the first-order position-momentum:
        DO J=HI2,II2
         RR2 = X_CELL(NODE,CELL,J) + XP(IP,J) * RR1
         X_CELL(NODE,CELL,J) = RR2
        END DO
        IF (COUNT(CELL) .EQ. II3) IFIRST(CELL) = I
        I = I + HI2
        IF (I .LE. N) CHANGING = T(I,HI2) .GT. T(I-HI2,HI2)
        NEXTNODE = CHANGING .OR. I .GT. N
        IF (NEXTNODE) CALL NONDEGENERACY
       END DO
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
 
 
 
 
      SUBROUTINE LINKTREE
       INTEGER CELL, J, L, IP
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
       INTEGER II5, II4, II3, II2, II1
       INTEGER*2 HI2, HI1
       PARAMETER (II5 = 4, HI2 = 2, II4 = 0, II3 = 3, II2 = 1)
       PARAMETER (II1 = 8, HI1 = 1)
       REAL RR2, RR1
       LOGICAL LL1
       REAL RR13 (86), RR12 (86), RR11 (86), RR10 (86), RR9 (86), RR8 (
     X   86), RR7 (86), RR6 (86), RR5 (86), RR4 (86), RR3 (86)
 
       DO CELL=HI1,II1
 
 
 
        NP(NODE,CELL) = COUNT(CELL)
!        make centroid:
        LL1 = COUNT(CELL) .GE. II2
        IF (LL1) THEN
         RR1 = II2 / CELLMASS(NODE,CELL)
!          resume cell centroide:
         DO J=HI1,II3
          RR2 = X_CELL(NODE,CELL,J) * RR1
          X_CELL(NODE,CELL,J) = RR2
         END DO
        END IF
       END DO
 
       DO CELL=HI1,II1
!
        IF (COUNT(CELL) .GT. II2) THEN
         DOWN(NODE,CELL) = NEXT
         AREST(NEXT) = 0.5 * AREST(NODE)
         NEXT = NEXT + HI1
        ELSE
         DOWN(NODE,CELL) = II4
         IF (COUNT(CELL) .EQ. II2) THEN
          L = IFIRST(CELL)
          IP = T(L,HI2)
          S(IP) = AREST(NODE)
          T(L,II5) = II2
          LABEL(NODE,CELL) = IP
         END IF
        END IF
       END DO
       NODE = NODE + HI1
       IF (NODE .GT. 40959) THEN
        PRINT *, 'LinkTree: ERROR: tree overflow'
        PRINT *, 'node=', NODE, '  nnode=', 40959
        STOP 
       END IF
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
 
 
 
 
      SUBROUTINE NONDEGENERACY
       INTEGER CELL, I, I1, IDEGEN
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
       COMMON /MOMENTA/ VMOMENT, RR1, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, RR2, U_OLD, U_DOT, RR3, U_DOT_VISC
     X   , LAMBDA, RR4, VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, RR5, KSEED, LSEED, RR6, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, RR7, Q_ACC, RR8, P_ACC
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
       COMMON /TREE_INTG/ NP, RR9, LABEL, RR10, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, RR11, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
       COMMON /NON_SINGULARITY/ ISINGULARITY
       COMMON /INDICES/ I, I1
C*$* Padding ( RR11, RR10, RR9, RR8, RR7, RR6, RR5, RR4, RR3, RR2, RR1 )
       INTEGER*2 HI2, HI1
       INTEGER II2, II1
       PARAMETER (HI2 = 0, II2 = 8, HI1 = 1, II1 = 0)
       REAL RR11 (86), RR10 (86), RR9 (86), RR8 (86), RR7 (86), RR6 (86)
     X   , RR5 (86), RR4 (86), RR3 (86), RR2 (86), RR1 (86)
       IDEGEN = II1
       DO CELL=HI1,II2
        IF (COUNT(CELL) .EQ. HI2) THEN
         IDEGEN = IDEGEN + HI1
        END IF
       END DO
       IF (FLOAT (IDEGEN) .EQ. 7.) THEN
        CALL CLR
        ISINGULARITY = ISINGULARITY + HI1
        CALL FINDOCTANTS (I1,I - HI1)
        AREST(NODE) = 0.5 * AREST(NODE)
        I = I1
       ELSE
        ISINGULARITY = II1
       END IF
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
 
 
 
 
 
 
      SUBROUTINE FINDOCTANTS ( IA, IB )
       INTEGER CELL, I, IA, IB, J, IP
 
 
 
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
       COMMON /MOMENTA/ VMOMENT, RR1, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, RR2, U_OLD, U_DOT, RR3, U_DOT_VISC
     X   , LAMBDA, RR4, VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, RR5, KSEED, LSEED, RR6, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, RR7, Q_ACC, RR8, P_ACC
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
       COMMON /TREE_INTG/ NP, RR9, LABEL, RR10, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, RR11, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
       COMMON /NON_SINGULARITY/ ISINGULARITY
C*$* Padding ( RR11, RR10, RR9, RR8, RR7, RR6, RR5, RR4, RR3, RR2, RR1 )
       INTEGER*2 HI3, HI2, HI1
       INTEGER II3, II2
       PARAMETER (HI3 = 3, II3 = 3, HI2 = 1, HI1 = 2, II2 = 0)
       INTEGER II1
       REAL RR11 (86), RR10 (86), RR9 (86), RR8 (86), RR7 (86), RR6 (86)
     X   , RR5 (86), RR4 (86), RR3 (86), RR2 (86), RR1 (86)
       DO I=IA,IB
        CELL = II2
        II1 = T(I,HI1)
        DO J=HI2,II3
         CELL = CELL * HI1
         BOX(II1,J) = BOX(II1,J) * HI1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
         IF (BOX(II1,J) .GT. 1.) THEN
          BOX(II1,J) = BOX(II1,J) - HI2
          CELL = CELL + HI2
         END IF
 
         IF (BOX(II1,J) .GT. 1.) THEN
          PRINT *, 'FindOctants: error'
          PRINT *, 'box=', BOX(II1,J)
          STOP 
         END IF
 
        END DO
        T(I,HI3) = CELL
       END DO
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
 
 
 
 
      SUBROUTINE CLR
       INTEGER CELL, J
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
       COMMON /MOMENTA/ VMOMENT, RR1, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, RR2, U_OLD, U_DOT, RR3, U_DOT_VISC
     X   , LAMBDA, RR4, VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, RR5, KSEED, LSEED, RR6, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, RR7, Q_ACC, RR8, P_ACC
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
       COMMON /TREE_INTG/ NP, RR9, LABEL, RR10, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, RR11, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
C*$* Padding ( RR11, RR10, RR9, RR8, RR7, RR6, RR5, RR4, RR3, RR2, RR1 )
       INTEGER II3, II2, II1
       INTEGER*2 HI1
       PARAMETER (II3 = 3, II2 = 0, II1 = 8, HI1 = 1)
       REAL RR11 (86), RR10 (86), RR9 (86), RR8 (86), RR7 (86), RR6 (86)
     X   , RR5 (86), RR4 (86), RR3 (86), RR2 (86), RR1 (86)
 
       DO CELL=HI1,II1
 
 
 
        COUNT(CELL) = II2
        CELLMASS(NODE,CELL) = 0.
       END DO
 
       DO CELL=HI1,II1
        DO J=HI1,II3
         X_CELL(NODE,CELL,J) = 0.
        END DO
       END DO
       NEXTNODE = .FALSE.
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
!
!
      SUBROUTINE NOLEAVES
       INTEGER K, I, J
       LOGICAL NOLEAF
!
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
       COMMON /MOMENTA/ VMOMENT, RR1, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, RR2, U_OLD, U_DOT, RR3, U_DOT_VISC
     X   , LAMBDA, RR4, VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, RR5, KSEED, LSEED, RR6, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, RR7, Q_ACC, RR8, P_ACC
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
       COMMON /TREE_INTG/ NP, RR9, LABEL, RR10, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, RR11, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
C*$* Padding ( RR11, RR10, RR9, RR8, RR7, RR6, RR5, RR4, RR3, RR2, RR1 )
       INTEGER*2 HI4, HI3, HI2, HI1
       INTEGER II4, II3, II2
       PARAMETER (HI4 = 2, HI3 = 3, HI2 = 0, II4 = 4, II3 = 1)
       PARAMETER (HI1 = 1, II2 = 0)
       EXTERNAL MPPTID, MPTEPA, MPTXPA, MPPFKD, MPOFRK, PKNOLEAVES_
       INTEGER MPPTID, MPPFKD, II16(0:9)
       CHARACTER AA3*26
       INTEGER II15, II1
       REAL RR11 (86), RR10 (86), RR9 (86), RR8 (86), RR7 (86), RR6 (86)
     X   , RR5 (86), RR4 (86), RR3 (86), RR2 (86), RR1 (86)
       AUTOMATIC II15
       DATA II16(0)/1/ 
       DATA II16(1)/0/ 
       DATA II16(2)/26/ 
       DATA AA3/';SPH_tree.F;NOLEAVES;388;;'/ 
       EQUIVALENCE (AA3, II16(3))
       INTEGER II17(0:9)
       CHARACTER AA4*26
       EXTERNAL PKNOLEAVES_1
       DATA II17(0)/0/ 
       DATA II17(1)/0/ 
       DATA II17(2)/26/ 
       DATA AA4/';SPH_tree.F;NOLEAVES;408;;'/ 
       EQUIVALENCE (AA4, II17(3))
       INTEGER II20, II19, II18
       PARAMETER (II20 = 334, II19 = 0, II18 = 200)
       II15 = MPPTID ()
!
!      copia somente as linhas nao-folha
!      e dispensa a coluna dos "leaf-flags"(4):
!
       K = II2
       DO I=HI1,N
        IF (T(I,II4) .NE. II3) THEN
         K = K + HI1
         DO J=HI1,II4
          T(K,J) = T(I,J)
         END DO
        END IF
       END DO
       N = K
       IF (K .EQ. HI2) RETURN 
!
!      calcula o novo parametro de peano promovendo os octantes para o 
!      status de no':
!
       IF (K .GT. II18 .AND. MPPFKD () .EQ. II19) THEN
        CALL MPOFRK (PKNOLEAVES_,II16,K)
       ELSE
C!!!!! PARALLEL REGION IF (K .GT. 200) SHARED (K,T) LOCAL (II1,I)
        CALL MPTEPA (II15)
C!!!!! PARALLEL DO 
        DO I=HI1,K
 
         II1 = T(I,HI1) * 8. + T(I,HI3)
         T(I,HI1) = II1
 
 
 
        END DO
C!!!!! END PARALLEL DO NOWAIT
        CALL MPTXPA (II15)
C!!!!! END PARALLEL REGION 
       END IF
!
!      ordena o array t por heap-sort dispensando a coluna dos 
!      octantes (col. 3):
!
       CALL HSORT
!
!      reduz o parametro de peano a fim de evitar estouro de inteiro:
!
       K = II2
       DO I=HI4,N
        IF (T(I,HI1) .NE. T(I-HI1,HI1)) K = K + HI1
        T(I-HI1,HI3) = K
       END DO
       IF (N .GT. II20 .AND. MPPFKD () .EQ. II19) THEN
        CALL MPOFRK (PKNOLEAVES_1,II17)
       ELSE
C!!!!! PARALLEL REGION IF (N .GT. 334) SHARED (N,T) LOCAL (I)
        CALL MPTEPA (II15)
C!!!!! PARALLEL DO 
        DO I=HI4,N
         T(I,HI1) = T(I-HI1,HI3)
        END DO
C!!!!! END PARALLEL DO NOWAIT
        CALL MPTXPA (II15)
C!!!!! END PARALLEL REGION 
       END IF
       T(HI1,HI1) = II2
!
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
      SUBROUTINE PKNOLEAVES_1 ( MPPID, MPPNPR )
       AUTOMATIC II10
       INTEGER II10
       INTEGER MPPID
       INTEGER MPPNPR
       INTEGER*2 HI12
       PARAMETER (HI12 = 1)
       INTEGER FOUR2
       PARAMETER (FOUR2 = 4)
       AUTOMATIC II14
       INTEGER II14
       INTEGER*2 HI32
       PARAMETER (HI32 = 3)
       INTEGER NN2
       PARAMETER (NN2 = 40960)
       AUTOMATIC II13
       INTEGER II13
       INTEGER*2 HI41
       PARAMETER (HI41 = 2)
       INTEGER N1
       AUTOMATIC II12
       INTEGER II12
       AUTOMATIC I2
       INTEGER I2
       INTEGER T4(40960,4)
       COMMON /TOTAL_NUMBER/ N1
       COMMON /T/ T4
C!!!!! PARALLEL REGION IF (N .GT. 334) SHARED (N,T) LOCAL (I)
       INTEGER II21
       PARAMETER (II21 = 1)
       II13 = N1 - HI41 + II21
       II14 = (II13 + MPPNPR - II21) / MPPNPR
       II10 = HI41 + MPPID * II14
       II12 = MIN (N1, II10 + (II14 - II21))
C!!!!! PARALLEL DO 
       DO I2=II10,II12
        T4(I2,HI12) = T4(I2-HI12,HI32)
       END DO
C!!!!! END PARALLEL DO NOWAIT
C!!!!! END PARALLEL REGION 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
!
!      calcula o novo parametro de peano promovendo os octantes para o 
!      status de no':
!
      SUBROUTINE PKNOLEAVES_ ( MPPID, MPPNPR, K1 )
       AUTOMATIC II11
       INTEGER II11
       INTEGER*2 HI31
       PARAMETER (HI31 = 3)
       EXTERNAL MPPTID
       INTEGER MPPTID
       AUTOMATIC II5
       INTEGER II5
       AUTOMATIC II8
       INTEGER II8
       INTEGER FOUR1
       PARAMETER (FOUR1 = 4)
       INTEGER*8 MPPFOD1(0:7,0:127)
       INTEGER*2 HI11
       PARAMETER (HI11 = 1)
       EXTERNAL MPPIOA
       INTEGER MPPIOA
       INTEGER MPPID
       AUTOMATIC II7
       INTEGER II7
       INTEGER NN1
       PARAMETER (NN1 = 40960)
       AUTOMATIC I1
       INTEGER I1
       INTEGER MPPNPR
       AUTOMATIC II9
       INTEGER II9
       AUTOMATIC II6
       INTEGER II6
       INTEGER K1
       INTEGER MPPFOA1
       COMMON /MPPFOA/ MPPFOA1
       INTEGER T3(40960,4)
       COMMON /MPPFOD/ MPPFOD1
       COMMON /T/ T3
       INTEGER II23, II22
       PARAMETER (II23 = 1, II22 = 0)
       II9 = MPPTID ()
       IF (MPPFOA1 .LT. II22) MPPFOD1(MPPIOA (II23),II9) = %LOC (K1)
C!!!!! PARALLEL REGION IF (K .GT. 200) SHARED (K,T) LOCAL (II1,I)
       II7 = K1 - HI11 + II23
       II8 = (II7 + MPPNPR - II23) / MPPNPR
       II5 = HI11 + MPPID * II8
       II6 = MIN (K1, II5 + (II8 - II23))
C!!!!! PARALLEL DO 
       DO I1=II5,II6
 
        II11 = T3(I1,HI11) * 8. + T3(I1,HI31)
        T3(I1,HI11) = II11
 
 
 
       END DO
C!!!!! END PARALLEL DO NOWAIT
C!!!!! END PARALLEL REGION 
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
!
!
      SUBROUTINE HSORT
!
!            method used: heapsort (num. rec. 1986)
!
       INTEGER L, IR, K, I, J
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
       COMMON /MOMENTA/ VMOMENT, RR1, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, RR2, U_OLD, U_DOT, RR3, U_DOT_VISC
     X   , LAMBDA, RR4, VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, RR5, KSEED, LSEED, RR6, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, RR7, Q_ACC, RR8, P_ACC
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
       COMMON /TREE_INTG/ NP, RR9, LABEL, RR10, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, RR11, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
C*$* Padding ( RR11, RR10, RR9, RR8, RR7, RR6, RR5, RR4, RR3, RR2, RR1 )
       INTEGER II4, II3
       INTEGER*2 HI2, HI1
       PARAMETER (II4 = 2, II3 = 1, HI2 = 2, HI1 = 1)
       INTEGER II2, II1
       REAL RR11 (86), RR10 (86), RR9 (86), RR8 (86), RR7 (86), RR6 (86)
     X   , RR5 (86), RR4 (86), RR3 (86), RR2 (86), RR1 (86)
!
!
       L = N / HI2 + HI1
       IR = N
   10  CONTINUE
       IF (L .GT. II3) THEN
        L = L - HI1
        DO K=HI1,HI2
         ROW(K) = T(L,K)
         II1 = T(HI1,K)
         T(L,K) = II1
        END DO
       ELSE
        DO K=HI1,HI2
         ROW(K) = T(IR,K)
         II2 = T(HI1,K)
         T(IR,K) = II2
        END DO
        IR = IR - HI1
        IF (IR .EQ. II3) THEN
         DO K=HI1,HI2
          T(HI1,K) = ROW(K)
         END DO
         RETURN 
        END IF
       END IF
       I = L
       J = L * II4
   20  IF (J .LE. IR) THEN
        IF (J .LT. IR) THEN
         IF (T(J,HI1) .LT. T(J+HI1,HI1)) J = J + HI1
        END IF
        IF (ROW(HI1) .LT. T(J,HI1)) THEN
         DO K=HI1,HI2
          T(I,K) = T(J,K)
         END DO
         I = J
         J = J * II4
        ELSE
         J = IR + HI1
        END IF
        GO TO 20
       END IF
       DO K=HI1,HI2
        T(I,K) = ROW(K)
       END DO
       GO TO 10
      END
C     KAP/Digital_UA_F      4.0 k3011126 980529 o5r3so3  17-Jul-2000 18:47:44
!
!
      SUBROUTINE FORBID
!
       INTEGER CELL
!
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
       COMMON /MOMENTA/ VMOMENT, RR1, ANGMOM
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
       COMMON /THERMAL_ENERG/ U, U_O, RR2, U_OLD, U_DOT, RR3, U_DOT_VISC
     X   , LAMBDA, RR4, VMONAG
       DIMENSION NSEED(40960), MSEED(40960), LSEED(40960), KSEED(40960),
     X    JSEED(40960), ISEED(40960)
 
 
 
 
       COMMON /SEEDS/ ISEED, JSEED, RR5, KSEED, LSEED, RR6, MSEED, NSEED
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
 
 
       COMMON /FORCES/ G, RR7, Q_ACC, RR8, P_ACC
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
       COMMON /TREE_INTG/ NP, RR9, LABEL, RR10, DOWN
       DIMENSION GNODE(327672,3)
 
       COMMON /GNODES/ GNODE
       DIMENSION XP(40960,3)
 
 
       COMMON /POSITIONS/ XP
       DIMENSION V(40960,3)
       COMMON /VELOCITIES/ V
       DIMENSION V_OLD(40960,3), V_O(40960,3)
 
       COMMON /V0/ V_O, RR11, V_OLD
       DIMENSION ACCEL(40960,3)
       COMMON /RESULT_ACC/ ACCEL
       DIMENSION BOX(40960,3)
 
 
 
 
 
 
 
 
       COMMON /BOXY_COORDINATES/ BOX
 
C*$* Padding ( RR11, RR10, RR9, RR8, RR7, RR6, RR5, RR4, RR3, RR2, RR1 )
       INTEGER*2 HI2, HI1
       INTEGER II1
       PARAMETER (HI2 = 0, II1 = 8, HI1 = 1)
       REAL RR11 (86), RR10 (86), RR9 (86), RR8 (86), RR7 (86), RR6 (86)
     X   , RR5 (86), RR4 (86), RR3 (86), RR2 (86), RR1 (86)
!
 
       DO CELL=HI1,II1
 
 
 
        NEXT = DOWN(NODE,CELL)
        IF (NEXT .GT. HI2) PERMIT(NEXT) = .FALSE.
       END DO
      END
!
!
# 536 "SPH_tree.F"
 
 
!
!
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
