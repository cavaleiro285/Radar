      SUBROUTINE RYLM(D_COLAT,D_LON,I_ORDER,D_YLMVAL)
c
c     $Revision: 1.2 $
c     The initial version of this subroutine was written by RADEX, INC.
c     for use on a VAX/VMS system.  
c     
c     Revisions for use with UNIX systems were written by KBB at 
c     The Johns Hopkins Univ. Applied Physics Laboratory.
c
c     These subsequent revisions have been maintained using the 
c     Revision Control System (RCS), and a log of all the changes
c     will be found at the end of the comments.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C     Purpose:
C
C         This subroutine computes the array of real spherical harmonic
C         function values for a given colatitude and longitude up to order
C         I_ORDER
C
C
C     Input Arguments:
C
C         D_COLAT             - Double Precision  - The colatitude of the
C                               point for which the spherical harmonic
C                               Y(L,M) will be calculated
C
C         D_LON               - Double Precision  - The longitude of the
C                               point for which the spherical harmonic
C                               Y(L,M) will be calculated
C
C         I_ORDER             - Integer  - The order of the spherical
C                               harmonic function expansion. The total
C                               number of terms computed will be
C                               (I_ORDER + 1) * (I_ORDER + 1)
C
C     Output Argument:
C
C         D_YLMVAL            - Double Precision array of spherical harmonic
C                               functions at the point (D_COLAT, D_LON)
C
C
C     Local Variables:
C
C     Constants:  None
C
C     Subroutines Required: None
C
C     Files Used at Compile Time:  None
C
C     Files Used at Run Time:  None
C
C     Databases Accessed:  None
C
C     Warnings:  None
C
C     Revision History:
C
C     Written by Radex, Inc., 3 Preston Court, Bedford, MA 01730
C
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     RCS   Revision History
c
c     $Log:	rylm.f,v $
c Revision 1.2  94/10/14  10:50:45  10:50:45  baker (Kile Baker S1G)
c Added the SAVE instruction to make the variables static.
c 
c Revision 1.1  94/10/12  15:24:21  15:24:21  baker (Kile Baker S1G)
c Initial revision
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SAVE   !Make all variables static

      DOUBLE PRECISION D_YLMVAL(1)
      DOUBLE PRECISION D_COLAT
      DOUBLE PRECISION D_LON
C
      COMPLEX*16 Q_FAC
      COMPLEX*16 Q_VAL
C
C     LOCAL VARIABLES
C
      DOUBLE PRECISION COS_THETA
      DOUBLE PRECISION SIN_THETA
      DOUBLE PRECISION COS_LON
      DOUBLE PRECISION SIN_LON
C
      COS_THETA =  DCOS(D_COLAT)
      SIN_THETA =  DSIN(D_COLAT)
C
      COS_LON   =  DCOS(D_LON)
      SIN_LON   =  DSIN(D_LON)
C
      Q_FAC     = -SIN_THETA * DCMPLX(COS_LON,SIN_LON)
C
C     GENERATE ZONAL HARMONICS (Y(L,0), L = 1, I_ORDER)
C     USING RECURSION RELATIONS
C
      D_YLMVAL(1) = 1.0D0
      D_YLMVAL(3) = COS_THETA
C
      DO 50 L = 1, I_ORDER - 1
      L0 = (L - 1) * L + 1
      L1 = L * (L + 1) + 1
      L2 = (L + 1) * (L + 2) + 1
C
      C1 = (2.0D0 * L + 1.0D0)/(L + 1.0D0)
      C2 = DFLOAT(L)/(L + 1.0D0)
C
      D_YLMVAL(L2) = C1 * COS_THETA * D_YLMVAL(L1) - C2 * D_YLMVAL(L0)
C
50    CONTINUE
C
C     GENERATE Y(L,L) FOR L = 1 TO I_ORDER
C
      Q_VAL = Q_FAC
C
      D_YLMVAL(4) =  DREAL(Q_VAL)
      D_YLMVAL(2) = -DIMAG(Q_VAL)

      DO 100 L = 2, I_ORDER
C
      Q_VAL = (2.0D0 * L - 1.0D0) * Q_FAC * Q_VAL
C
      L0 = L * L
      L1 = L0 + 2 * L + 1
      L2 = L0 + 1
C
      D_YLMVAL(L1) =  DREAL(Q_VAL)
      D_YLMVAL(L2) = -DIMAG(Q_VAL)
C
100   CONTINUE
C
C     GENERATE Y(L+1,L) TERMS
C
      DO 150 L = 2, I_ORDER
C
      LA  = L  * L
      L1  = LA + 2 * L
C
      LB  = LA - 2 * (L - 1)
      L2  = L1 - 2 * (L - 1)
C
      FAC = 2.0D0 * L - 1.0D0
C
      D_YLMVAL(L1) = FAC * COS_THETA * D_YLMVAL(LA)
      D_YLMVAL(L2) = FAC * COS_THETA * D_YLMVAL(LB)
C
150   CONTINUE
C
C     GENERATE Y(L+1,L), YL(L+2,L) ..., YL(I_ORDER,L)
C
      DO 200 M = 1, I_ORDER - 2
      LX = M + 2
C
      LA0 =  (LX - 1)**2
      LA1 =   LX * LX - 1
      LA2 =  (LX + 1)**2 - 2
C
      LB0 =   LA0 - 2 * M
      LB1 =   LA1 - 2 * M
      LB2 =   LA2 - 2 * M
C
      DO 180 L = LX, I_ORDER
C
C
      C1  = DFLOAT(2*L-1)/DFLOAT(L-M)
      C2  = DFLOAT(L+M-1)/DFLOAT(L-M)
C
      D_YLMVAL(LA2) = C1 * COS_THETA * D_YLMVAL(LA1) -
     X C2 * D_YLMVAL(LA0)
C
      D_YLMVAL(LB2) = C1 * COS_THETA * D_YLMVAL(LB1) -
     X C2 * D_YLMVAL(LB0)
C
      LA0 = LA1
      LA1 = LA2
      LA2 = LA2 + 2 * L + 2
C
      LB0 =   LA0 - 2 * M
      LB1 =   LA1 - 2 * M
      LB2 =   LA2 - 2 * M
C
180   CONTINUE
C
200   CONTINUE
500   CONTINUE
      RETURN
      END
