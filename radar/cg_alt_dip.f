      SUBROUTINE CG_ALT_DIP(ALTITUDE,ALATITUDE,IFLAG,R_LATITUDE,I_ERR)
C
C     $Revision: 1.3 $
C
c     The initial version of this program was written by Radex, Inc. 
c     on 11 Oct. 1994.
c
c     The original version was written for a VAX/VMS system.
c
c     Revisions for UNIX systems have been made by KBB at 
c     The Johns Hopkins Univ. Applied Physics Laboratory.
c
c     These additional revisions have been made utilizing the 
C     Revision Control System (RCS).
c     
c     The log of all the revisions made via RCS will be found at the end
c     of the comments:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C     Except for low altitudes, the corrected geomagnetic latitude
C     is discontinuous in the vicinity of the magnetic equator. This
C     continuity would result in poor convergence of a spherical
C     harmonic expansion. To avoid this problem, it was necessary to
C     perform the spherical harmonic fit and expansion with respect
C     to dipole coordinates at the table altitudes, and to compute
C     an appropriate latitude adjustment. This routine computes the
C     required adjustment.
C
C     COMPUTES LATITUDE ADJUSTMENT FOR ALTITUDE FOR USE WITH
C     SPHERICAL HARMONIC EXPANSION
C
C     Input Arguments:
C
C         ALTITUDE    - Single Precision - Input height in km above
C                       earth mean radius
C
C         ALATITUDE   - Single Precision - Altitude Corrected Dipole
C                       Latitude [IFLAG = 1] or Corrected Geomagnetic
C                       Latitude [IFLAG = 2]
C
C         IFLAG       - Integer - Direction Flag for Computation
C
C                       = 1 convert altitude corrected dipole latitude
C                           to corrected geomagnetic coordinates
C                       = 2 convert corrected geomagnetic coordinates
C                           to altitude corrected dipole coordinates
C
C     Output Arguments:
C
C         R_LATITUDE  - Single Precision - Corrected Geomagnetic Latitude
C                       (IFLAG = 1) or Altitude Corrected Dipole Coordinates
C                       (IFLAG = 2)
C
C         I_ERR       - Integer - Error Flag
C
C                       = 0 normal return
C 
C                       = 1 indicates that the input coordinates for
C                         the IFLAG = 2 option are invalid for the
C                         input altitude.  
C                         
C     Constants:
C
C         ERADIUS     - Single Precision - Earth Radius in km
C
C         EPS         - Single Precision - limit parameters, to avoid
C         UNIM        - Single Precision - computational singularities
C
C     Internal Variables:
C
C         RA, R0      - Single Precision - for intermediate results
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     $Log:	cg_alt_dip.f,v $
c Revision 1.3  94/10/17  12:37:07  12:37:07  baker (Kile Baker S1G)
c added error code return value.
c 
c Revision 1.2  94/10/14  10:37:06  10:37:06  baker (Kile Baker S1G)
c This version adds the SAVE instruction to make the variables static.
c 
c Revision 1.1  94/10/12  15:19:40  15:19:40  baker (Kile Baker S1G)
c Initial revision
c 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SAVE   !make all variables static

      DATA ERADIUS / 6371.2 /
      DATA EPS     / 1.0E-9 /
      DATA UNIM    / 0.9999999 /
C
      I_ERR = 0
C
      IF (IFLAG .EQ. 1) THEN
C
      RA = (COS(0.017453293 * ALATITUDE))**2
      IF (RA .LT. EPS) RA = EPS
      R0 = (1.0 + ALTITUDE/ERADIUS)/RA
      IF (R0 .LT. UNIM) R0 = UNIM
C
      R_LATITUDE = 57.295779713 * SIGN(ACOS(SQRT(1.0/R0)), ALATITUDE)
      RETURN
C
      ELSE
C
      RA = (1.0 + ALTITUDE/ERADIUS) * (COS(0.017453293 * ALATITUDE))**2
      IF (RA .GT. UNIM) THEN 
          RA = UNIM
          I_ERR = 1
      ENDIF
      R_LATITUDE = 57.295779713 * SIGN(ACOS(SQRT(RA)), ALATITUDE)
C
      ENDIF
      RETURN
      END
