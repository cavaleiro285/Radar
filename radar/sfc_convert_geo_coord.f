      SUBROUTINE SFC$$CONVERT_GEO_COORD(R_LAT_IN,
     $                                  R_LON_IN,
     $                                  R_HEIGHT_IN,
     $                                  R_LAT_OUT,
     $                                  R_LON_OUT,
     $                                  I_FLAG,
     $                                  I_ERROR)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     $Revision: 1.3 $
c    The initial version of this subroutine was written by RADEX, INC.
c     for use on VAX/VMS systems.
c
c     Subsequent revisions for UNIX systems have been made by KBB at
c     The Johns Hopkins Univ. Applied Physics Laboratory.  These revisions
c     have been managed using the Revision Control System (RCS) and a
c     log of the revisions will be found at the end of the comments.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C     Purpose:
C
C         This subroutine uses a set of spherical harmonic coefficients to
C         perform a coordinate conversion between geographic and corrected
C         geomagnetic coordinates (or vice versa).
C
C
C
C
C
C
C         The spherical harmonics for the geographic to corrected geomagnetic
C         coordinate system correspond to a conversion from geographic to
C         the corrected geomagnetic coordinates (expressed in terms of
C         centered dipole coordinates at the input altitude), and are then
C         transformed to ground level.
C
C         The spherical harmonic coefficients used for the inverse are
C         computed relative to centered dipole coordinates at the input
C         altitude. The input CGM coordinates are converted to equivalent
C         values in this coordinate sytem using the inverse altitude
C         algorithm before the evaluation of the geographic spherical
C         harmonic expansion.
C
C     Method:
C
C         This subroutine uses a five-step process in converting a position
C         in one coordinate system to a position in the other.
C         The five steps are as follows:
C
C         1.  The appropriate spherical harmonic coefficients for the
C             coordinate conversion are computed for the input altitude.
C
C         2.  The appropriate coordinates for use in the spherical harmonic
C             expansion are computed. For the geographic ==> corrected
C             geomagnetic coordinate version, these are geographic colatitude
C             and longitude. For the inverse coordinate conversion, the
C             input coordinates are first converted into equivalent dipole
C             coordinates at the input altitude.
C
C         3.  The Cartesian coordinates of the unit vector in the desired
C             coordinate system are then computed using the appropriate
C             spherical harmonic expansion.
C
C         4.  For the geographic ==> corrected geomagnetic coordinate
C             conversion, the dipole-equivalent coordinates at input
C             altitude are converted to their cgm cartesian coordinates.
C
C         5.  Standard trigonometric identities are used to compute the
C             latitude and longitude in the desired output coordinate system.
C
C     Input Arguments:
C
C         R_LAT_IN     - REAL*4 - The input latitude in degrees.
C                        This could be either geographic latitude
C                        or corrected geomagnetic latitude.  The
C                        acceptable range of values is (-90 to +90
C                        degrees).
C
C         R_LON_IN     - REAL*4 - The input longitude in degrees east.
C                        This could be either geographic longitude
C                        or corrected geomagnetic longitude.  The
C                        acceptable range of values is (0 to 360 degrees).
C
C         R_HEIGHT_IN  - REAL*4 - The height in kilometers for which the
C                        coordinate transformation will be accomplished.
C                        The acceptable range of values is (0 km to
C                        2000 km)
C
C         I_FLAG       - INTEGER - The flag that indicates which way
C                        the conversion will proceed.
C
C                        = 1  will convert geographic to corrected
C                                 geomagnetic coordinates
C
C                        = 2  will convert corrected geomagnetic
C                                 to geographic coordinates
C
C     Output Arguments:
C
C         R_LAT_OUT   - REAL*4 - The output latitude in degrees.
C                       This could be either the geographic latitude
C                       or corrected geomagnetic latitude.  This
C                       value will be between -90 and +90 degrees.
C
C         R_LON_OUT   - REAL*4 - The output longitude in degrees.
C                       This could be either geographic longitude or
C                       corrected geomagnetic longitude.  This
C                       value will be between 0 and 360 degrees.
C
C         I_ERROR     - INTEGER - The error flag
C
C                       =  0  normal processing
C                       = -2  R_HEIGHT_IN is outside the allowable range
C                       = -4  I_FLAG value is invalid
C                       = -8  R_LAT_IN is outside the allowable range
C                       = -16 R_LON_IN is outside the allowable range
C                       = -32 Magnitude of the "unit vector" of the
C                             target coordinate system deviates
C                             significantly (+/- 10% or more) from 1.
C                       = -64 For altitudes > 0, for the corrected 
C                             geomagnetic to geographic coordinate
C                             conversion there are regions for which
C                             for low latitudes there are regions in which
C
C     Local Variables:
C
C         D_CINT(,,)       - Double Precision Array (3-D) -
C                            Contains the spherical harmonic coefficients
C                            interpolated to the input height, R_HEIGHT_IN.
C
C         D_COEF(,,,)      - Double Precision Array (3-D) -
C                            Contains the spherical harmonic coefficients
C                            used to compute the Cartesian coordinates of
C                            unit vector in the target coordinate system.
C
C                            First index: sp. harm. coeff. index
C                            Second       x, y, z components of unit vector
C                            Third        altitude indices 0, 300, 1200 km
C                            Fourth       direction of conversion index
C
C                            coefficients for a given altitude have the form
C
C                            a0 + h a1 + h * h a2 where h = alt/1000 [km]
C
C         D_COLAT_INPUT    - Double Precision - colatitude (radians) in the
C                            input coordinate system
C
C         D_COLAT_OUTPUT   - Double Precision - colatitude (radians) in the
C                            output coordinate system
C
C         D_LON_INPUT      - Double Precision - longitude (radians) in the
C                            input  input system
C
C         D_LON_OUTPUT     - Double Precision - longitude (radians) in the
C                            output coordinate system
C
C         D_R              - Double Precision - magnitude of the
C                            unit radius vector in the target coordinate
C                            system.  The target coordinate system is
C                            the system to which this subroutine is
C                            converting the latitude and ongitude.
C
C         D_X              - Double Precision - the X-component of the
C                            unit radius vector in the target coordinate
C                            system
C
C         D_Y              - Double Precision - the Y-component of the
C                            unit radius vector in the target coordinate
C                            system
C
C         D_Z              - Double Precision - the Z-component of the
C                            unit radius vector in the target coordinate
C                            system
C
C         D_YLMVAL         - Double Precision - the array of spherical
C                            harmonic basis functions evaluated at
C                            a particular colatitude and longitude.
C
C         K                - Integer  - first index of the D_CINT interpolated
C                            spherical harmonic coefficient array
C
C         L                - Integer  - order of each spherical harmonic
C                            coefficient and spherical harmonic function
C
C         M                - Integer  - zonal index of the spherical
C                            harmonic functions
C
C         R_HEIGHT_GRID()  - Real*4 Array (1-D)  - Contains the heights
C                            of the four grid levels in kilometers
C
C         R_HEIGHT_OLD()   - Real*4 Array (1-D)  - Variable containing
C                            previous height (km) used to determine whether
C                            the interpolation to compute D_CINT() from
C                            D_COEF() needs to be done.
C
C     Constants:
C
C         I_MAX_LENGTH     - Integer  - the length of the D_COEF(,,,)
C                            and D_YLMVAL spherical harmonic function
C                            arrays, corresponding to the order of the
C                            spherical harmonic expansion used.
C                            I_MAX_LENGTH = (I_ORDER + 1) * (I_ORDER + 1)
C
C         I_NUM_AXES       - Integer  - the number of axes in a Cartesian
C                            coordinate system (3)
C
C         I_NUM_FLAG       - Integer  - the number of coordinate trans-
C                            formation flags available (2)
C
C         I_NUM_LEVEL      - Integer  - the number of grid levels
C                            available.  This is the number of
C                            distinct heights for which there are
C                            spherical harmonic coefficients available
C
C         I_ORDER          - Integer  - the order of the spherical harmonic
C                            expansion used.
C
C         DEGRAD           - Double Precision  - the multiplicative
C                            conversion factor for converting degrees
C                            to radians
C
C         PI               - Double Precision  - the well-known
C                            mathematical constant
C
C     Subroutines Required:
C
C         CINTRP           - This subroutine interpolates the Y(L,M)
C                            coefficients for the specified height
C                            R_HEIGHT_IN using the coefficients corresponding
C                            to the heights in the R_HEIGHT_GRID array.
C
C         RYLM              - This subroutine returns the spherical harmonic
C                            function value array for a given colatitude and
C                            longitude.  The spherical harmonics are returned
C                            in the D_YLMVAL array.
C
C
C     Files Used at Compile Time:  None
C
C
C     Files Used at Run Time:  None
C
C
C     Databases Accessed:  None
C
C
C     Warnings:
C
C
C
C
C
C
C
C
C
C
C
C
C     Revision History:
C
C             Original Version was converted to a standard library
C             subroutine for SFC$$ by MFS on  14 OCT 1992.
C
C             Revisions by Radex, Inc. (Contractor to USAF Phillips Lab.,
C             Geophysics Dir.) Bedford, MA 01730, 11 Oct 94
C             Tel: 617-275-6767
C
C        The Following Revisions have been made by Radex, Inc. 10/94
C
C        (1) The set of spherical harmonic function values are computed 
C            recursively from a single call to subroutine RYLM rather 
C            than using a function call. Computationally this method 
C            is much faster than computing them as needed.
C
C        (2) The new spherical harmonic coefficients are based upon 
C            the IGRF 90 magnetic field model. The order of the
C            spherical harmonic expansion used here is 10 (previously 
C            was = 4).
C
C        (3) The coefficients represent a quadratic fit to the spherical
C            harmonic fits at 0, 300 and 1200 km altitude. The allowable  
C            range of altitudes is 0 - 2000 km. The previous version used
C            polynomial interpolations of the coefficients at 0, 150,
C            300 and 450 km. 
C  
C        (4) The spherical harmonic coefficients were computed in the
C            applicable target coordinate system, resulting in a
C            considerable simplification of the code. The need to perform
C            multiple coordinate system rotations has been eliminated,
C            and the unnecessary code removed.
C
C
C
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     RCS  Revision Log
c
c     $Log:	sfc_convert_geo_coord.f,v $
c Revision 1.3  94/10/17  12:35:32  12:35:32  baker (Kile Baker S1G)
c added error code -64 to indicate invalid magnetic coordinates specified
c as input.  This also requires a change in the call to cg_alt_dip
c 
c Revision 1.2  94/10/14  10:53:36  10:53:36  baker (Kile Baker S1G)
c Added the SAVE instruction to make variables static
c 
c Revision 1.1  94/10/12  15:28:38  15:28:38  baker (Kile Baker S1G)
c Initial revision
c 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C         Subroutine begins:
C
      SAVE    !make all variables static

      PARAMETER (PI     = 3.141592656589793D0 )
      PARAMETER (DEGRAD = 1.745329251994330D-2)
C
C         Set the order used in this algorithm to 10.
C
      PARAMETER (I_ORDER = 10)
C
C         Parameterize the maximum values of the indices of the
C         spherical harmonic coefficient array D_COEF( , , , ) and the
C         interpolated spherical harmonic coefficient array, D_CINT( , , ).
C
      PARAMETER (I_NUM_TERMS = 121)
      PARAMETER (I_NUM_AXES  =   3)
      PARAMETER (I_NUM_LEVEL =   3)
      PARAMETER (I_NUM_FLAG  =   2)
C
      DOUBLE PRECISION D_CINT(I_NUM_TERMS,I_NUM_AXES,I_NUM_FLAG)
      DOUBLE PRECISION D_COEF(I_NUM_TERMS,I_NUM_AXES,I_NUM_LEVEL,
     $                        I_NUM_FLAG)
C
      DOUBLE PRECISION D_COLAT_INPUT
      DOUBLE PRECISION D_COLAT_OUTPUT
      DOUBLE PRECISION D_LON_INPUT
      DOUBLE PRECISION D_LON_OUTPUT
C
      DOUBLE PRECISION D_X
      DOUBLE PRECISION D_Y
      DOUBLE PRECISION D_Z
      DOUBLE PRECISION D_R
C
      DOUBLE PRECISION D_YLMVAL(I_NUM_TERMS)
C
      REAL R_HEIGHT_OLD(2)
      REAL R_HEIGHT_GRID(3)
C
C         Initialize R_HEIGHT_OLD() to impossible values.
C
      DATA R_HEIGHT_OLD(1)/-1./
      DATA R_HEIGHT_OLD(2)/-1./
C
C         Data in the heights of each grid level (in kilometers).
C
      DATA R_HEIGHT_GRID(1)/    0.0 /
      DATA R_HEIGHT_GRID(2)/  300.0 /
      DATA R_HEIGHT_GRID(3)/ 1200.0 /
C
C         The long list of data statements that follows assigns values
C         to the spherical harmonic coefficients used to generate the
C         Cartesian coordinates of the unit vector in the target coordinate
C         system.  The target coordinate system is the system to which
C         the algorithm is converting and is the coordinate system
C         corresponding to the output variables, R_LAT_OUT and R_LON_OUT.
C
C         The coefficient set is stored in a 4-dimensional array with
C         indices defined as follows:
C
C         First Index  - Represents the number of terms used in the
C                        spherical harmonic expansion. Equal to
C                        (I_ORDER + 1) * (I_ORDER + 1)
C
C         Second Index - Represents the Cartesian coordinate a particular
C                        coefficient will be used to generate.  Indices are
C                        defined as follows:
C
C                        1  - X-coordinate (the X-axis points from the
C                             center of the earth to where the prime
C                             meridian crosses the equator.
C
C                        2  - Y-coordinate (the Y-axis points from the
C                             center of the earth to 90 degrees east
C                             of the X-axis)
C
C                        3  - Z-coordinate (the Z-axis points from the
C                             center of the earth to the north pole)
C
C         Third Index  - Represents the terms of a quadratic fit to
C                        altitude (independent variable h =
C                        altitude [km]/1200) for a given spharical
C                        harmonic coefficient. The indices are defined
C                        as follows:
C
C                        1  - Constant term: corresponds to the fit at
C                             0 km altitude
C                        2  - linear term
C                        3  - quadratic term
C
C         Fourth Index - Represents the direction of the coordinate
C                        transformation.  Indices are defined as follows:
C
C                        1  - Conversion of geographic coordinates to
C                             corrected geomagnetic coordinates.
C
C                        2  - Conversion of corrected geomagnetic coordinates
C                             to geographioc coordinates.
C
C
C
      DATA (D_COEF(J,1,1,1),J=1,62)/     2.4117280D-02, -9.1550411D-01,
     X  -1.9019433D-01, -3.2418227D-01, -1.4506193D-03,  9.3644609D-03,
     X   4.9722584D-02, -9.1896281D-03,  8.2716243D-04,  2.4630113D-05,
     X  -2.1346006D-04,  1.0353331D-02, -1.2298884D-02, -3.2030517D-03,
     X   2.3756025D-03, -3.7396049D-05,  1.1331450D-05,  1.5909722D-04,
     X  -3.8503450D-04,  2.0339766D-03, -7.2406229D-03,  6.7275292D-04,
     X  -5.1827021D-05, -2.2394954D-04,  1.6745100D-05,  7.8966323D-07,
     X  -1.0273538D-05,  2.1841898D-05,  2.9150452D-04,  1.4269359D-04,
     X   4.9633155D-03,  1.6216637D-03,  6.4826128D-05,  2.7644659D-05,
     X   1.6867073D-07, -1.3499715D-06, -1.8185712D-08, -1.4685169D-08,
     X   2.6087997D-06,  4.7465067D-06, -8.3635064D-05,  3.8744218D-04,
     X  -1.5170154D-04,  2.0840154D-04, -1.2029954D-04,  7.9573035D-08,
     X  -2.0502275D-07, -7.4134185D-08,  3.7222804D-08,  1.4878152D-09,
     X   1.5977409D-10, -1.4793917D-07,  5.7438694D-07,  7.0402877D-06,
     X  -5.9281545D-05, -5.0258326D-05,  1.0312565D-03, -3.5443916D-04,
     X  -6.9048556D-06,  6.6304603D-06,  6.6853249D-07, -4.8745925D-08/
C
      DATA (D_COEF(J,1,1,1),J=63,121)/  -1.3211130D-08, -2.2807147D-09,
     X   8.2060775D-11, -8.2801234D-10,  4.2498656D-09, -2.2510827D-08,
     X  -1.1286149D-07,  3.0842367D-06,  8.3453806D-06, -1.1673996D-04,
     X  -1.2366731D-03, -1.8475305D-04,  1.2308199D-05,  1.7611476D-06,
     X   4.5774641D-08, -3.4938537D-08, -2.4594692D-09, -2.1257687D-10,
     X   1.8226687D-11, -3.5538728D-12, -7.8068685D-12, -2.1289482D-11,
     X   1.7802465D-09,  1.4187161D-08, -1.6902279D-07, -1.3659845D-06,
     X   1.5938850D-05,  8.6180091D-05, -1.3966055D-04,  4.2085540D-05,
     X   1.6334565D-05,  9.6060624D-08, -1.2897418D-07, -3.0749935D-09,
     X  -2.9973672D-10,  9.5597378D-11, -2.8009276D-12, -1.7638741D-12,
     X   2.4913307D-13, -2.6456595D-13, -7.4749691D-12, -2.7323756D-11,
     X  -1.4551349D-10,  9.0110217D-09,  4.4544009D-08, -8.3670893D-07,
     X  -7.6162246D-06,  5.6352683D-05, -4.2708908D-04,  1.7212712D-04,
     X   1.9861296D-06, -9.4803132D-07, -4.6494190D-08,  2.7863011D-09,
     X   5.7809422D-10,  6.3650214D-11, -8.4172228D-13, -1.9063722D-12,
     X  -9.9255412D-14/
C
      DATA (D_COEF(J,2,1,1),J=1,62)/     1.3298518D-02,  3.1393449D-01,
     X   2.2524234D-03, -9.4924268D-01,  3.6391845D-03, -6.7630690D-03,
     X   2.8790576D-02, -3.9837571D-02, -2.4569589D-03, -5.5342960D-04,
     X   1.6891812D-03, -1.3765321D-03, -4.2246270D-02, -5.8498630D-04,
     X   4.5372898D-03,  4.5145620D-04, -3.2186035D-06, -2.4828775D-04,
     X  -6.3907504D-04, -2.9950592D-04,  1.0102714D-02,  3.8544848D-03,
     X  -1.3948746D-03,  1.1526202D-04, -4.1784105D-05,  2.0645877D-06,
     X  -6.5648014D-06,  7.6564923D-06, -1.3269768D-04,  3.1939340D-05,
     X   6.6774892D-03, -7.9203687D-04, -1.2584176D-04,  4.0374993D-05,
     X  -1.6106501D-05,  1.0782209D-06, -1.8299012D-07,  7.7510068D-07,
     X  -6.5532978D-07,  1.5203247D-06, -8.1713400D-05,  1.6116168D-04,
     X   2.3544710D-03,  1.0548488D-03, -1.4381759D-04, -4.1479719D-06,
     X   6.7090583D-06,  1.2423827D-06,  8.5119467D-09,  2.8999223D-09,
     X  -7.0148301D-09, -1.4563629D-07,  1.6188493D-07,  2.4652583D-06,
     X  -1.2425323D-05, -1.3305083D-04,  1.7872781D-03, -4.6427340D-05,
     X  -9.0547901D-05,  2.4886915D-06,  1.2297799D-06, -9.2664417D-08/
C
      DATA (D_COEF(J,2,1,1),J=63,121)/  -8.1191266D-09, -5.6673741D-10,
     X  -2.3988639D-11,  3.9317567D-10,  1.1940247D-08, -1.7319798D-08,
     X  -1.6795664D-07,  6.1427408D-07,  4.9107629D-06, -5.0111305D-05,
     X  -2.8675231D-04, -5.4514651D-04,  9.7877759D-06,  7.8660559D-06,
     X   1.0753474D-07, -8.5811201D-08, -4.1105915D-09, -2.6404285D-11,
     X   1.6312478D-10, -1.6987162D-12, -2.5267027D-11, -1.0077528D-10,
     X   1.7395792D-09,  1.0418091D-08, -7.6546736D-08, -1.5606082D-06,
     X   1.7934281D-06,  8.4982461D-05, -2.7028213D-03,  5.2960122D-06,
     X   4.8154512D-05,  4.1228279D-07, -3.8408559D-07, -1.1434520D-08,
     X   2.2677269D-09,  3.0124617D-10,  1.5738918D-12, -7.9283666D-12,
     X   3.3453118D-13,  1.0741274D-13, -5.7312057D-12, -9.3468937D-11,
     X  -2.7696021D-10,  9.2416036D-09,  8.1999228D-08, -2.4523024D-07,
     X  -1.1636874D-05,  3.5803428D-05, -4.6611632D-04,  4.5218960D-04,
     X   7.9812854D-06, -3.1846437D-06, -8.7964902D-08,  1.0411850D-08,
     X   1.0714832D-09, -1.5077260D-12, -1.2461785D-11, -1.0747204D-12,
     X   1.1536701D-13/
C
      DATA (D_COEF(J,3,1,1),J=1,62)/    -2.6744177D-02, -1.6495709D-01,
     X   9.8545363D-01,  9.5272707D-03,  1.2882348D-03,  2.1436797D-02,
     X   4.7925086D-02,  1.6444218D-02, -2.3173807D-02,  7.8641517D-04,
     X  -2.0105738D-03, -2.1913168D-03, -5.3007384D-03, -1.6855191D-02,
     X  -3.1442814D-03,  2.5637341D-03, -2.0996507D-04,  1.7954882D-04,
     X  -2.4684249D-05, -1.4424547D-03,  3.8256031D-03, -2.7253843D-03,
     X  -1.0266528D-03,  5.4690652D-05, -1.0578892D-05, -1.5071085D-06,
     X  -1.7385603D-05,  2.1967887D-05,  5.8006368D-04,  6.8680001D-04,
     X  -1.1427431D-02, -3.4697930D-03,  5.5035386D-04,  7.1009792D-05,
     X  -4.9921255D-06, -8.7413628D-06,  2.9106898D-07,  5.7070927D-07,
     X  -5.3721820D-07, -5.4842860D-05,  3.8177190D-05,  5.8330711D-04,
     X  -1.7601356D-02,  3.4516925D-03,  6.0764930D-04, -6.3433357D-05,
     X  -1.1003497D-05,  3.4878323D-07,  8.6640427D-07, -1.7506029D-08,
     X  -7.6394167D-08,  3.0965466D-07,  3.8343660D-06, -5.6047398D-06,
     X  -1.4607244D-04, -4.2356918D-05,  5.9269441D-03,  3.8265209D-03,
     X  -1.7113011D-04, -4.7195451D-05,  3.0205700D-07,  4.8547688D-07/
C
      DATA (D_COEF(J,3,1,1),J=63,121)/   4.7185968D-08, -9.1545335D-10,
     X   2.6379366D-10,  1.6421919D-09, -2.3409636D-08, -1.2683903D-07,
     X   5.5461202D-07,  1.2007828D-05, -3.8561956D-05, -4.5297286D-04,
     X   1.3670295D-02, -7.2123927D-04, -3.2064395D-04,  2.6434355D-06,
     X   2.7633439D-06,  1.6636927D-07, -1.1817240D-08, -5.3735771D-09,
     X  -4.2739139D-10, -6.9863429D-12,  6.9422060D-11,  9.0415670D-10,
     X   1.7990525D-09, -7.3416693D-08, -6.7986010D-07,  3.0403821D-06,
     X   5.5038899D-05,  1.9834253D-04, -3.8023991D-03, -2.2887223D-03,
     X   2.3687620D-06,  1.7854100D-05,  5.6204941D-07, -6.6539170D-08,
     X  -1.1552980D-08, -7.1054407D-10,  2.1434920D-10,  3.8197601D-11,
     X  -1.5952804D-12, -4.4426175D-12, -3.7185207D-11,  2.3163910D-10,
     X   4.0581788D-09,  1.5035872D-08, -1.3149858D-07, -4.2061534D-06,
     X   1.2832394D-05,  2.3055051D-04, -9.0962295D-03, -1.4435260D-04,
     X   1.3345211D-04,  2.8380929D-06, -6.1320968D-07, -5.4317021D-08,
     X  -1.1633049D-09,  5.1416814D-10,  6.4636229D-11, -2.0826507D-12,
     X  -1.3003818D-12/
C
      DATA (D_COEF(J,1,2,1),J=1,62)/    -4.7289327D-03, -5.7258676D-03,
     X   1.6615938D-03,  1.5081695D-03,  6.7743997D-04, -1.7499320D-03,
     X  -1.0406227D-02,  8.3523813D-04,  4.6034614D-04, -3.1406013D-05,
     X   1.3454734D-04, -3.5019532D-03,  3.2213875D-03,  1.8008864D-03,
     X  -8.5541502D-04,  1.2934680D-05, -7.3175801D-06, -5.3462346D-05,
     X  -4.9738704D-05, -1.7505607D-03,  2.1726320D-03, -3.5224192D-05,
     X   1.5571984D-04,  1.3235667D-04, -9.3639637D-06,  1.4147928D-07,
     X   4.4545432D-06,  2.1042164D-08, -1.4186386D-04, -2.6190897D-04,
     X  -2.5257558D-03, -7.4225970D-04, -1.2741494D-04, -2.4849346D-05,
     X   1.5458895D-08,  1.1298206D-06, -5.8814348D-08,  2.0850622D-07,
     X  -1.7676013D-06,  1.7452853D-06,  1.7085042D-05, -5.2870244D-04,
     X  -1.8486941D-04, -3.8645448D-04,  8.5113865D-05,  3.1467192D-06,
     X   1.0239559D-06,  1.6403429D-08, -1.9582992D-08,  2.7726815D-09,
     X  -7.2634907D-09,  4.9903718D-08, -5.7300041D-07, -2.8679116D-06,
     X   6.6916079D-05,  1.7439939D-05, -1.6636089D-03,  9.9702956D-05,
     X  -7.4865848D-06, -2.4747365D-06, -1.6830927D-07,  2.8604086D-08/
C
      DATA (D_COEF(J,1,2,1),J=63,121)/   6.7792764D-09, -9.0600526D-10,
     X  -1.9674379D-10,  1.8517649D-09,  2.9574019D-09, -4.8417835D-09,
     X  -2.7576214D-07, -2.7831362D-06,  2.4157882D-06,  9.0027087D-05,
     X   1.3629283D-04,  2.0265027D-04,  1.5624966D-05,  7.0905702D-07,
     X  -2.5238623D-07, -2.8267004D-09,  9.8123246D-11,  3.9792220D-10,
     X   2.4110692D-10,  1.1077871D-12, -3.9751961D-11, -2.0554488D-10,
     X   1.0957159D-09,  1.4211372D-08,  8.1704019D-08, -4.7645386D-07,
     X  -1.1579079D-05,  7.1373896D-06, -9.8462447D-04,  9.1835565D-05,
     X  -1.2042405D-05, -1.0312396D-06, -3.3797729D-08, -3.8533988D-09,
     X   2.0912914D-09,  2.8985376D-10,  3.4507828D-11, -9.4177783D-12,
     X   1.4099860D-13,  1.8507756D-12,  5.5389774D-12, -1.7261056D-10,
     X  -1.0911620D-09,  6.2602201D-09,  9.9301013D-08,  1.2485971D-07,
     X  -6.6784685D-06, -6.0237714D-05,  8.8384609D-04,  7.6237998D-06,
     X  -1.2948942D-05, -9.8772840D-07,  1.2411452D-08,  7.1748699D-09,
     X   6.4361913D-10, -9.3241715D-11, -2.0710694D-11,  1.3231263D-12,
     X   5.9142169D-13/
C
      DATA (D_COEF(J,2,2,1),J=1,62)/    -5.5799008D-04,  3.6561273D-03,
     X   3.2566986D-05, -3.1409610D-03, -7.7318619D-04,  1.9035953D-03,
     X  -1.0073237D-02,  7.7013033D-03,  2.0236279D-03,  1.8199637D-04,
     X  -7.5713726D-04,  4.8234571D-04,  1.4143968D-02,  1.6800170D-03,
     X  -1.5392473D-03, -3.6772958D-04,  5.7247715D-06,  1.4770192D-04,
     X   1.3296988D-04,  6.5786644D-04, -4.1384378D-03, -8.6704295D-04,
     X   4.9250152D-04, -6.8556152D-05,  2.8012319D-05, -1.6785345D-06,
     X   5.5337045D-06,  7.5444285D-07,  1.8337101D-05, -3.5007667D-04,
     X  -2.0364638D-03,  5.4266700D-04,  9.9105915D-06, -2.1841195D-05,
     X   1.1788365D-05, -2.7434388D-07,  1.5865390D-07, -7.6605432D-07,
     X  -1.4829012D-07,  5.2384478D-06,  6.5114895D-05,  6.0622917D-05,
     X  -1.5147778D-03, -1.0038580D-03,  8.6345014D-05,  1.1808655D-05,
     X  -4.2130242D-06, -1.0763017D-06, -9.0542216D-08, -1.3341296D-09,
     X   1.9936686D-08,  1.5493958D-07, -2.6184486D-07, -2.9682377D-06,
     X   1.3875937D-05,  5.3022553D-05, -1.7772886D-03, -1.8997055D-04,
     X   7.5980829D-05,  2.7956265D-06, -7.4773967D-07, -2.6842451D-08/
C
      DATA (D_COEF(J,2,2,1),J=63,121)/  -4.7925765D-09,  5.6074786D-09,
     X  -1.7085824D-10, -4.5050701D-10, -8.9825852D-09,  1.1588413D-08,
     X  -3.7455608D-09, -1.5595490D-08, -6.2865323D-06,  1.0494890D-04,
     X  -9.6064078D-04,  3.8718127D-04,  2.4992890D-05, -5.5371249D-06,
     X  -4.8846070D-07,  2.4942791D-08,  1.2178350D-08,  9.5077802D-10,
     X  -3.8552019D-10,  2.0427213D-11,  1.9020576D-12, -1.9949512D-10,
     X  -1.3394737D-09, -2.2279013D-10,  7.8269624D-08,  4.0261442D-08,
     X   9.7858354D-06, -4.7089489D-05,  1.1005697D-03,  2.6856718D-04,
     X  -5.3304515D-06, -2.2412035D-06,  3.4882215D-08,  2.6644639D-08,
     X   1.4194562D-09, -4.5938647D-10, -5.8350485D-11,  9.8261422D-12,
     X  -1.2218350D-12,  2.8476854D-12,  1.5367601D-11,  5.9740424D-11,
     X  -6.1078745D-10, -4.7142564D-09,  2.7196015D-08, -1.4835232D-07,
     X  -4.2122330D-06,  2.0980651D-05,  1.5090662D-03,  1.9091546D-05,
     X  -2.3196023D-05, -4.5456254D-07,  9.7871664D-08,  1.7169559D-08,
     X  -1.3686575D-09, -1.6398283D-10,  8.9902520D-12,  2.2243604D-12,
     X   3.7167737D-13/
C
      DATA (D_COEF(J,3,2,1),J=1,62)/     4.9007943D-03, -2.3473771D-03,
     X  -8.2410452D-03, -2.0540049D-02, -1.2037815D-03, -4.2239204D-03,
     X  -1.3460616D-02,  7.5257013D-03,  2.8125723D-03,  3.9581561D-05,
     X   4.8595727D-04, -5.4909750D-04,  2.4881589D-02,  9.0184737D-03,
     X  -1.7921942D-03, -1.3057930D-03,  1.1707724D-04,  7.3247338D-05,
     X   3.8910885D-04, -7.3381651D-04,  9.5509031D-03, -1.0068509D-02,
     X  -4.2206727D-04,  3.7019125D-04,  4.7532394D-05, -9.7302083D-07,
     X  -1.2775280D-05, -6.6205281D-05,  1.7678798D-04,  1.4058138D-03,
     X  -2.9297936D-02, -5.9038812D-03,  1.3116393D-03,  2.0333285D-04,
     X  -2.6317749D-05,  3.6934446D-07,  3.2233584D-07,  1.7034614D-06,
     X  -3.6853554D-07, -6.5389431D-05, -2.1214320D-04,  1.2795180D-03,
     X  -2.1853132D-02,  7.3525132D-03,  1.0231971D-03, -9.2488516D-05,
     X  -1.6752993D-05, -2.3868900D-07, -5.3971517D-07, -6.4891606D-09,
     X   5.4731163D-09,  2.6576770D-07,  6.7834069D-06,  1.3071470D-05,
     X  -3.0865026D-04, -1.4763135D-03,  2.6963487D-02,  8.3540665D-03,
     X  -6.5707335D-04, -1.3230776D-04,  2.4873581D-06,  1.1826497D-06/
C
      DATA (D_COEF(J,3,2,1),J=63,121)/   8.8425774D-08, -4.9070865D-09,
     X   1.6376844D-09, -5.6262215D-09, -5.5721079D-08, -4.3712356D-07,
     X   4.6476196D-07,  5.0686189D-05,  9.5032180D-05, -1.4750655D-03,
     X   4.1263914D-02, -4.7138977D-03, -9.6111595D-04,  3.0100143D-05,
     X   8.7045375D-06,  3.5837453D-07, -3.2895985D-08, -7.5636885D-09,
     X  -1.1914713D-09,  3.9763543D-11,  5.6487938D-10,  1.9883008D-09,
     X   2.9649703D-09, -1.3192780D-07, -3.7602238D-06, -4.8375338D-06,
     X   2.6269402D-04,  1.2171072D-03, -1.5386785D-02, -9.2956299D-03,
     X   1.8474962D-04,  8.0170982D-05,  1.5166145D-06, -3.7263799D-07,
     X  -3.9384316D-08,  3.3229333D-11,  1.3695162D-10,  1.0008400D-11,
     X  -1.0636087D-14, -3.3507387D-11, -1.0381607D-10,  1.1576963D-09,
     X   1.5718734D-08,  1.4676558D-07, -1.4566911D-07, -3.1502491D-05,
     X  -3.5239943D-05,  1.5926798D-03, -4.9176424D-02,  9.6897690D-04,
     X   7.2542595D-04,  5.7601776D-06, -3.8539630D-06, -2.6740945D-07,
     X   4.0474164D-09,  2.2832232D-09,  1.4354259D-10,  8.5266984D-12,
     X   1.5453544D-12/
C
      DATA (D_COEF(J,1,3,1),J=1,62)/     6.3147483D-04,  1.4533142D-03,
     X  -4.4753384D-04,  5.3215312D-05, -2.6554034D-04,  3.0282524D-04,
     X   2.1266811D-03,  2.1928014D-04, -2.5825778D-04,  2.2125843D-05,
     X  -6.2919785D-05,  5.8102914D-04, -3.6808123D-04, -6.2314885D-04,
     X   1.5804347D-04, -4.2619877D-06,  3.5177837D-08,  1.5706689D-05,
     X   3.8283735D-05,  5.1290743D-04, -5.9964326D-04, -1.2427437D-04,
     X  -5.7540382D-05, -3.4651670D-05,  2.8111924D-06, -1.9657595D-07,
     X  -1.4170412D-06, -2.0599356D-06,  4.2100914D-05,  9.0514726D-05,
     X   6.4370705D-04,  1.7032037D-04,  6.6284384D-05,  9.8595058D-06,
     X   5.9433393D-08, -3.8627587D-07,  4.4800959D-08, -9.7712237D-08,
     X   5.2055078D-07, -1.5139280D-06,  1.5226408D-06,  2.2220207D-04,
     X  -1.7201318D-04,  1.7858027D-04, -2.6334870D-05, -1.5412795D-06,
     X  -5.9654331D-07, -1.0388982D-08,  4.6299240D-09, -2.7345684D-09,
     X   5.5502805D-09,  1.2922104D-08,  1.8615058D-07,  6.1889017D-08,
     X  -2.4355275D-05, -9.6786891D-06,  5.4245758D-04,  2.2643305D-05,
     X   7.7579790D-06,  3.0090350D-07, -8.2570435D-08, -1.2440160D-08/
C
      DATA (D_COEF(J,1,3,1),J=63,121)/  -6.9041605D-11,  1.4654161D-09,
     X   1.3431332D-10, -1.1097394D-09, -3.5400134D-09,  1.4711166D-08,
     X   2.0203117D-07,  6.6657275D-07, -8.1060736D-06, -4.9456094D-06,
     X   4.1871838D-06, -5.4131211D-05, -1.1420328D-05, -1.2887252D-06,
     X   1.1669210D-07,  9.8164296D-09,  1.2030311D-09, -1.6219581D-10,
     X  -2.1501794D-10,  1.7010725D-12,  3.3847584D-11,  1.2875068D-10,
     X  -1.4935914D-09, -1.4200355D-08,  1.2717960D-08,  7.5767088D-07,
     X   2.3986065D-06, -4.9661023D-05,  8.7604878D-04, -6.9500547D-05,
     X   1.3912654D-06,  6.0092664D-07,  6.8005128D-08,  6.0972137D-09,
     X  -1.2684113D-09, -2.5686267D-10, -2.8935964D-11,  9.7441826D-12,
     X  -3.3034305D-13, -1.2320219D-12,  6.8483232D-13,  1.3941372D-10,
     X   6.4747734D-10, -8.6470276D-09, -6.5438356D-08,  4.1533067D-07,
     X   6.3368912D-06,  3.0843603D-05, -6.2727981D-04, -7.4551855D-05,
     X   6.4722570D-06,  1.0428546D-06,  2.2196034D-08, -4.7169808D-09,
     X  -6.9595415D-10,  2.7983036D-11,  1.4327862D-11,  1.0407694D-13,
     X  -4.2646012D-13/
C
      DATA (D_COEF(J,2,3,1),J=1,62)/    -7.1924769D-04, -7.6505303D-04,
     X  -4.0411988D-04,  2.5656303D-03,  4.6019850D-05, -4.1282927D-04,
     X   3.2386256D-03, -8.2514350D-04, -9.2921779D-04, -1.7294206D-05,
     X   1.7649553D-04, -9.0616851D-05, -2.3290482D-03, -1.0256849D-03,
     X   2.2020744D-04,  1.4786264D-04, -4.0702893D-06, -4.0274417D-05,
     X  -1.7666397D-05, -2.1458737D-04,  5.2417726D-04, -1.3855530D-04,
     X  -7.9648058D-05,  3.3400158D-05, -1.0140995D-05,  6.8696150D-07,
     X  -2.0533720D-06, -8.5012100D-07,  1.3274094D-06,  1.6206194D-04,
     X  -1.5915224D-04, -9.3370784D-05,  3.1655528D-05,  5.6519516D-06,
     X  -4.4478958D-06, -5.0096809D-08, -5.4323533D-08,  2.9128644D-07,
     X   1.0627677D-07, -3.1334019D-06, -3.0360173D-05, -5.9081736D-05,
     X   5.3890410D-04,  4.8102418D-04, -2.8654895D-05, -6.4879320D-06,
     X   1.1423041D-06,  3.5699890D-07,  6.5858686D-08, -6.8223599D-10,
     X  -9.4882364D-09, -5.8098637D-08,  7.7950557D-08,  1.3551147D-06,
     X  -1.6286299D-06,  3.1270809D-05,  7.7534520D-04,  9.5252940D-05,
     X  -2.5939401D-05, -2.2902682D-06,  1.3323096D-07,  3.8534752D-08/
C
      DATA (D_COEF(J,2,3,1),J=63,121)/   7.5314422D-09, -3.7081904D-09,
     X   1.5868338D-10,  6.8865413D-11,  2.2906318D-09,  3.3831830D-09,
     X   6.1720177D-08, -4.5259187D-07, -1.1145087D-06, -3.7130451D-05,
     X   6.6104777D-04, -6.3726351D-05, -1.6570835D-05,  1.0229566D-06,
     X   2.7419310D-07,  1.3658855D-08, -6.9391927D-09, -6.5705578D-10,
     X   2.0079191D-10, -1.4406867D-11,  1.4548828D-11,  1.9299422D-10,
     X  -7.6769124D-11, -4.6973626D-09,  1.1161216D-08,  3.9276514D-07,
     X  -8.2908514D-06,  2.6268359D-07,  1.3748821D-04, -1.4930451D-04,
     X  -1.2722931D-05,  1.2042251D-06,  1.1697883D-07, -1.1547547D-08,
     X  -2.0283063D-09,  2.0458179D-10,  4.0435241D-11, -3.2931387D-12,
     X   8.2801247D-13, -2.1465745D-12, -8.6091465D-12,  5.2831972D-12,
     X   4.8798016D-10, -2.4308104D-09, -4.0476042D-08,  5.1525565D-07,
     X   4.1080994D-06, -4.7833222D-05, -6.2524593D-04, -1.6313847D-04,
     X   9.3026043D-06,  1.3690529D-06, -1.6149690D-08, -1.6790127D-08,
     X   4.0306267D-10,  1.2311640D-10, -6.4705687D-13, -1.0231761D-12,
     X  -3.2959162D-13/
C
      DATA (D_COEF(J,3,3,1),J=1,62)/    -5.1525554D-05, -4.5326068D-04,
     X   6.4050773D-03,  3.2599446D-03,  7.9041659D-04, -5.7202917D-05,
     X   4.1028406D-03, -8.4106349D-03, -2.8113208D-04, -1.4978150D-04,
     X   3.2049828D-04,  6.9600659D-04, -1.9459912D-02, -4.8466571D-03,
     X   2.0777083D-03,  5.1491131D-04, -3.1814251D-05, -1.3250456D-04,
     X  -2.0491141D-04,  1.2530431D-03, -1.1814976D-02,  8.8456535D-03,
     X   8.9373438D-04, -2.9764986D-04, -4.3751564D-05,  2.3219904D-06,
     X   1.7236264D-05,  2.5527703D-05, -3.7285138D-04, -9.6184024D-04,
     X   2.5493711D-02,  6.8445868D-03, -1.1273218D-03, -1.9833228D-04,
     X   1.7751574D-05,  3.2589906D-06, -4.2721087D-07, -1.4780645D-06,
     X   1.9815644D-06,  7.0771695D-05,  7.8454049D-05, -1.4188866D-03,
     X   2.6177303D-02, -6.3989868D-03, -1.1034652D-03,  7.8447432D-05,
     X   1.7446727D-05,  4.5567276D-07,  6.6100344D-08,  1.6338240D-08,
     X   1.4738721D-08, -5.2117834D-07, -6.4702351D-06, -9.5905498D-08,
     X   3.1490942D-04,  9.4851932D-04, -2.1799065D-02, -8.6840306D-03,
     X   4.7874354D-04,  1.2402282D-04, -3.9638327D-07, -1.0690800D-06/
C
      DATA (D_COEF(J,3,3,1),J=63,121)/  -1.1001798D-07,  3.6205674D-09,
     X  -1.5476749D-09,  4.9750447D-09,  6.0620407D-08,  3.2238152D-07,
     X  -1.0916626D-06, -4.3086736D-05, -1.9012181D-05,  1.4477220D-03,
     X  -3.9867046D-02,  3.1939096D-03,  9.1987596D-04, -1.3457671D-05,
     X  -7.7053163D-06, -4.4552014D-07,  2.3004516D-08,  9.0411689D-09,
     X   1.3285275D-09, -3.1447917D-11, -5.3980012D-10, -2.0673304D-09,
     X   1.9417040D-09,  1.8336417D-07,  2.9565207D-06, -1.3907584D-06,
     X  -2.3363164D-04, -9.0764438D-04,  9.4817681D-03,  8.5675669D-03,
     X  -5.9805991D-05, -6.9693741D-05, -2.1466315D-06,  2.8117856D-07,
     X   3.9050128D-08,  8.2693417D-10, -1.9508782D-10, -4.2384203D-11,
     X   1.8905382D-12,  3.0842158D-11,  9.5581632D-11, -1.3081734D-09,
     X  -1.5028881D-08, -8.9322457D-08,  3.5756447D-07,  2.5174304D-05,
     X  -9.3522596D-07, -1.4124326D-03,  4.3039258D-02,  1.0544147D-04,
     X  -6.3046533D-04, -1.2297937D-05,  3.0174012D-06,  2.6110424D-07,
     X   7.4051736D-10, -2.0126653D-09, -1.7220547D-10, -8.1059429D-12,
     X  -6.1834066D-13/
C
      DATA (D_COEF(J,1,1,2),J=1,62)/    -2.1073533D-02,  9.3707233D-01,
     X   5.8499372D-02, -3.0880580D-01, -1.9704620D-04, -3.5754359D-02,
     X  -4.0970361D-02,  6.0935837D-03, -1.2318387D-03, -1.2596043D-04,
     X   4.2361162D-03, -2.0520003D-03,  4.5659819D-02,  3.2891814D-03,
     X   3.5090457D-03, -5.6805858D-05, -3.2826485D-05, -2.6273204D-04,
     X  -1.3069066D-03,  7.6715352D-03, -1.5825763D-02, -3.1620401D-03,
     X  -1.3079137D-03, -2.4234815D-04, -5.5250839D-06,  3.2240538D-06,
     X   1.2449453D-05, -4.9066680D-05, -2.7381950D-04, -3.5867445D-03,
     X  -9.4200865D-03,  4.0936607D-04, -2.0189793D-04,  1.2328103D-04,
     X   7.0548163D-06,  4.2067568D-07, -1.7436014D-07, -8.3022882D-07,
     X  -1.4966540D-06,  1.6667738D-05,  7.8654719D-05,  1.4624707D-03,
     X   4.7699075D-03, -2.5415358D-05,  1.4414961D-05, -5.0177333D-06,
     X  -3.2508398D-06, -1.5198812D-06, -7.2132196D-08,  1.3837481D-08,
     X  -1.1316270D-08,  1.3084618D-07, -8.7147108D-07, -9.3700046D-06,
     X  -2.6906168D-05, -3.6585036D-04, -1.9257149D-03,  1.8344347D-04,
     X  -5.1837875D-05,  1.2566917D-05, -8.7817216D-08,  2.6038632D-07/
C
      DATA (D_COEF(J,1,1,2),J=63,121)/   2.8870652D-08,  1.0186694D-08,
     X  -2.8991315D-10, -1.1806052D-09,  6.2350229D-09, -6.4431150D-09,
     X   1.4153690D-06, -9.6892359D-07,  9.2026261D-05,  1.6639344D-04,
     X   6.1883321D-03, -8.8031298D-05,  1.1374018D-04, -1.2333881D-07,
     X   2.2776932D-07,  1.6815064D-08, -1.6630776D-08, -7.4513787D-10,
     X  -3.1040292D-10,  1.0759048D-11,  1.8233413D-11,  3.3137839D-10,
     X  -3.8665510D-09, -1.8696009D-08, -3.8241619D-07,  2.5360789D-06,
     X  -2.6839246D-05,  7.3906719D-04, -1.0196646D-03, -1.9842029D-04,
     X  -2.5789164D-05, -4.8462977D-06,  8.5426321D-08, -3.0484300D-08,
     X   4.7194047D-09,  3.1536566D-10,  6.2006382D-11,  1.2566258D-11,
     X  -8.9479052D-13,  9.5211125D-13, -1.0299902D-11,  1.2481553D-10,
     X   8.7269327D-10,  1.3401489D-08, -4.7293059D-08, -9.4675202D-07,
     X  -2.4691294D-05, -2.9376363D-04, -1.9604369D-03,  1.3122829D-04,
     X  -2.6718223D-05,  3.5895074D-06, -7.9910420D-09,  1.3054750D-08,
     X  -5.6195067D-10, -3.8280793D-11, -1.3745846D-11, -4.5351304D-12,
     X  -3.1062899D-13/
C
      DATA (D_COEF(J,2,1,2),J=1,62)/     1.2969429D-02,  3.2185529D-01,
     X  -1.7807459D-01,  9.4811321D-01,  3.2735940D-03,  5.8548353D-04,
     X   3.4814138D-02,  2.1875887D-02,  1.7260920D-03, -2.4887957D-04,
     X   1.0367934D-03,  7.6866425D-03, -1.0449116D-02,  6.2805458D-03,
     X  -8.2490035D-04, -6.2694611D-04, -1.3846131D-05, -3.8690834D-04,
     X   1.9605294D-04, -7.8933678D-04, -1.0458941D-02,  1.5969877D-03,
     X  -2.1953647D-04,  4.4028319D-05,  6.1684506D-05,  2.6062128D-06,
     X   2.6400852D-06, -3.8056840D-05, -3.1816375D-04,  4.1835787D-04,
     X   1.2368666D-03,  9.2801394D-04,  3.7682479D-04, -2.8927484D-05,
     X   1.1957947D-05, -4.1747636D-06, -2.5331546D-07, -1.1012750D-07,
     X  -2.6501033D-06,  6.9518047D-06,  9.8281300D-05, -3.1678642D-04,
     X  -4.0076466D-03,  2.5948965D-04, -1.2781845D-04,  5.7961759D-06,
     X  -4.1806911D-06, -6.6837965D-07,  7.3348852D-08,  1.1803946D-08,
     X  -3.9818915D-09, -4.4487505D-08,  7.0327476D-07, -1.1564068D-05,
     X   8.2027396D-05, -3.0167259D-04,  2.3237631D-03,  1.0663980D-04,
     X   7.8142843D-05, -6.9273927D-06, -8.1460083D-07, -5.2548288D-08/
C
      DATA (D_COEF(J,2,1,2),J=63,121)/  -2.0193562D-08,  4.0343592D-09,
     X  -1.7707958D-10,  5.5237439D-10,  5.9598173D-09,  2.1610394D-08,
     X  -2.1764173D-07,  1.4215424D-06, -8.8709604D-06,  6.4211774D-04,
     X  -2.3059752D-04, -2.6142370D-04,  6.4255728D-07, -5.1752692D-06,
     X   4.8702010D-07, -1.6348601D-08,  1.0080102D-09,  4.8465492D-10,
     X   9.1179204D-11, -2.3502954D-12,  3.3547237D-11,  5.3759744D-10,
     X   2.7549439D-09,  2.0771712D-08, -2.5541183D-07,  1.2773711D-06,
     X  -3.3786362D-05, -1.7081259D-04, -2.3402425D-03,  2.5935684D-04,
     X  -2.0263527D-05,  2.3855658D-06,  1.0472858D-07,  2.9717592D-09,
     X   2.7247595D-09,  1.6356405D-10, -6.8281086D-11,  1.5953548D-11,
     X  -4.9329343D-13, -4.2965836D-12,  6.7173277D-12, -2.0787142D-10,
     X   2.9343605D-10,  1.4459788D-08,  5.6263891D-08,  1.3282859D-06,
     X   1.1889973D-05, -1.6437678D-04,  4.4543278D-04,  1.0478739D-04,
     X   4.2144470D-06,  1.3570585D-06, -6.4764799D-08, -1.7975234D-08,
     X  -2.5873988D-10, -9.4966589D-11,  2.6663221D-11, -9.8294886D-13,
     X  -4.4817432D-13/
C
      DATA (D_COEF(J,3,1,2),J=1,62)/     2.9280647D-02,  6.0916195D-02,
     X   9.6944939D-01,  1.5454723D-01, -1.1907191D-02,  3.9126332D-03,
     X  -2.7334596D-02, -1.6513418D-02, -2.1185775D-02, -1.1594099D-03,
     X  -2.5254160D-03, -1.4283866D-02,  2.8520971D-03,  1.0820620D-02,
     X  -9.7270397D-04,  2.6372696D-03,  1.1940931D-04, -2.7505540D-05,
     X   1.7878529D-04,  1.4823643D-03,  4.7248564D-03,  6.0765988D-04,
     X  -3.8643147D-04,  5.3810191D-05, -1.9076522D-04, -1.7661581D-07,
     X   1.0965554D-05,  5.7671783D-06,  2.1528694D-04, -1.1052330D-03,
     X   8.6698549D-03,  1.1881680D-04,  1.4817271D-04, -2.6858807D-05,
     X  -2.7773825D-06,  4.2066880D-06,  5.9560391D-07, -3.0985353D-07,
     X   5.0967846D-06,  3.6596482D-05,  2.1910245D-04,  4.3742791D-03,
     X   7.8120244D-03, -1.1904319D-03,  2.4672637D-04, -7.9878423D-05,
     X   1.1520750D-06, -1.2607963D-06, -3.2780940D-07,  1.3529024D-09,
     X  -1.3277310D-08, -9.8313601D-08, -5.0885104D-06,  4.2202165D-08,
     X  -3.5499549D-04,  8.4696608D-04, -1.8714765D-02, -2.4728837D-04,
     X  -3.5252989D-04, -1.1915003D-05,  1.2190802D-06, -4.9794999D-08/
C
      DATA (D_COEF(J,3,1,2),J=63,121)/   1.0246780D-07,  7.2465866D-09,
     X  -3.5742411D-10,  1.9302590D-09,  4.9259578D-09,  1.5239631D-07,
     X   1.3265846D-07, -1.1927847D-05, -1.8809464D-05, -3.1102235D-03,
     X  -6.7745145D-04,  1.3345551D-03,  1.4387909D-06,  3.3270423D-05,
     X   1.8310377D-07,  1.1967528D-07, -1.6653840D-09, -3.7044113D-09,
     X  -5.7962796D-10, -1.6917086D-11, -1.7239005D-10, -7.4868265D-10,
     X  -2.6566471D-09, -1.5952345D-08,  1.3104028D-06,  2.2658367D-06,
     X   1.7244018D-04,  3.8595045D-04,  1.2076297D-02,  9.5914073D-08,
     X   1.4271874D-04, -3.9218121D-06, -4.6953922D-07, -4.0548212D-08,
     X  -6.8031807D-09,  1.7287059D-10,  1.5427921D-10,  2.7464678D-11,
     X   1.7356436D-12,  2.0330585D-12,  2.3412857D-11,  5.6474244D-11,
     X  -1.6153520D-09, -2.4845385D-08, -3.7584448D-07,  1.3126319D-06,
     X  -2.5089563D-05,  1.2924855D-03, -1.8402785D-03, -6.5564842D-04,
     X  -2.6834217D-05, -9.0882658D-06,  3.6513024D-08, -7.0120673D-09,
     X   3.0812308D-09,  2.7135038D-10,  1.3532902D-11, -1.0759806D-12,
     X  -7.1290514D-13/
C
      DATA (D_COEF(J,1,2,2),J=1,62)/     1.2806897D-03, -1.6653831D-03,
     X   2.1633202D-04, -1.7503586D-03,  8.2880881D-04,  4.7158277D-03,
     X   1.2438495D-02, -5.2127789D-04,  1.6331270D-03,  3.0186533D-04,
     X  -1.1783343D-03,  3.9241098D-03, -1.3161673D-02, -2.4412342D-03,
     X  -9.6484241D-04, -2.2421219D-04, -1.0514212D-05,  1.2079512D-04,
     X   8.9010250D-05, -1.0244536D-03,  2.9574197D-03,  3.7677812D-04,
     X   2.7795881D-04,  1.1393143D-04,  3.2776359D-06, -2.6717087D-06,
     X  -9.9194541D-06,  1.4512592D-05, -3.6777288D-05,  1.3850516D-04,
     X  -2.3583734D-03,  1.0426571D-03, -6.2654023D-05, -2.8224790D-05,
     X  -4.4155284D-06,  1.3148088D-06,  1.4622492D-07,  6.5325445D-07,
     X   5.6094623D-06, -2.0871296D-05,  9.8067561D-05, -3.2777557D-03,
     X   4.2410151D-03,  1.0842400D-03,  1.4042937D-04,  4.4454306D-05,
     X   4.9655797D-06,  8.7245553D-07,  7.8763546D-08, -1.2404060D-08,
     X   1.1969950D-08, -3.0362943D-07,  1.9173793D-06,  9.7392782D-06,
     X   2.2753348D-04,  1.2063351D-03,  1.2297933D-02, -5.9483935D-04,
     X   2.0344536D-04, -2.2654340D-05, -8.4066812D-07, -2.0264632D-07/
C
      DATA (D_COEF(J,1,2,2),J=63,121)/  -6.7306335D-09, -1.3421148D-08,
     X   1.8834073D-10,  4.2829595D-09, -7.4682486D-09,  3.7891842D-08,
     X  -7.5461247D-07,  6.3790125D-06, -1.7911107D-05,  2.0492743D-03,
     X   3.2451062D-04, -1.0929175D-03, -2.5932205D-05, -1.2843370D-05,
     X  -4.5194929D-07, -8.2906516D-08,  7.5373657D-09,  9.5079837D-11,
     X  -2.2564679D-10, -3.4253498D-11, -1.2071720D-11, -7.9098536D-10,
     X   8.6754377D-09,  3.4221272D-08,  6.2351690D-08, -1.3193914D-06,
     X  -2.0375167D-05,  1.1471410D-04, -3.3947382D-03, -1.9440659D-04,
     X   1.4528610D-05,  5.8533463D-06, -1.8149494D-08,  8.0027222D-08,
     X  -4.9492555D-09, -3.4435737D-11, -5.0834402D-11, -1.2598222D-12,
     X   2.7438894D-12, -5.7226502D-12,  1.2240768D-11, -4.6883450D-10,
     X  -6.7427250D-10,  1.0814811D-08, -2.7276450D-07,  1.1593640D-06,
     X  -6.8207133D-06,  2.5512366D-04, -2.9036102D-03,  5.5305006D-05,
     X  -3.6538886D-05, -2.8817829D-06, -5.0008823D-07,  5.5519630D-09,
     X  -8.5010707D-10, -2.5538438D-10,  6.9418479D-12,  1.3207220D-11,
     X   1.3966094D-12/
C
      DATA (D_COEF(J,2,2,2),J=1,62)/    -2.3213693D-03,  1.7746445D-03,
     X  -1.0928882D-03, -9.1149211D-03, -1.0353282D-03, -1.2255021D-03,
     X  -6.8703847D-03, -3.0012544D-03,  2.4052007D-05,  2.6016446D-04,
     X  -3.2802781D-04, -3.1447874D-03,  6.5718094D-03, -2.0550037D-03,
     X   7.6495340D-04,  4.5790092D-04,  1.1350697D-05,  1.3754741D-04,
     X   2.1880231D-04,  1.4570888D-03,  5.8146281D-03, -2.2451199D-03,
     X  -1.9821084D-04, -1.3686555D-04, -5.2853295D-05, -4.7779554D-06,
     X  -4.3567162D-06,  1.2746851D-05, -4.5632925D-05, -2.2838558D-04,
     X  -6.2469529D-03, -1.1214253D-03, -2.6958440D-04,  6.6037691D-06,
     X  -5.7768016D-06,  3.7846208D-06,  2.6275810D-07,  5.3042943D-07,
     X  -7.4047906D-07, -1.1578070D-06, -2.2216874D-04, -1.4351741D-03,
     X   1.5767535D-04,  4.4866462D-04,  1.9857856D-05,  7.2886253D-06,
     X   7.5467833D-06,  1.2310836D-07, -1.1459867D-08, -1.2626620D-08,
     X   1.8023274D-08,  2.5458794D-07,  9.8147839D-07,  9.5161106D-06,
     X   1.0861369D-04, -4.5577632D-05,  6.2083336D-03, -8.2946288D-05,
     X   2.9002751D-05,  1.7358801D-05,  7.3636224D-07,  2.5983848D-07/
C
      DATA (D_COEF(J,2,2,2),J=63,121)/   1.7705295D-08, -1.0624094D-08,
     X   4.7469907D-10,  9.3466694D-10, -2.6318534D-08, -6.1671844D-08,
     X   8.2277196D-07,  1.2125601D-05,  3.9482485D-05,  1.3923033D-03,
     X   2.1138687D-03, -3.8088228D-04,  3.7151236D-05, -8.8341997D-06,
     X  -1.3037450D-06, -4.3704581D-08, -8.1258452D-09, -1.2637780D-09,
     X  -8.9026312D-11, -1.1418516D-11, -9.2792680D-11, -2.1114858D-09,
     X   1.0725728D-09, -2.6215652D-08, -3.2094870D-07,  1.8220966D-06,
     X  -9.4999001D-05,  3.6564086D-04, -6.3185924D-03, -1.3383037D-04,
     X  -4.1362789D-05, -3.7922439D-06,  7.1678611D-07, -1.3173181D-08,
     X  -6.3228118D-09, -5.5188927D-10,  1.0286541D-10, -3.2546381D-11,
     X   1.5455146D-12,  8.9423884D-12,  5.5029384D-11,  4.2350646D-10,
     X  -1.6136814D-09,  2.4209531D-08, -3.9090395D-07, -9.3450175D-07,
     X   5.5985624D-05, -1.3731821D-03, -2.3338996D-03,  7.5110040D-04,
     X  -1.0185131D-04,  7.9090921D-06, -7.8049133D-07,  3.3604030D-09,
     X   3.0755081D-09, -5.1783056D-10, -7.4263789D-11,  9.5463905D-12,
     X   1.9794393D-12/
C
      DATA (D_COEF(J,3,2,2),J=1,62)/    -6.1805803D-03, -1.9339903D-02,
     X   1.0824039D-02,  9.9000243D-03,  1.3957692D-03,  9.4735791D-03,
     X   6.9229521D-03,  8.1049598D-04,  2.9129985D-03,  6.9161524D-04,
     X  -5.8682224D-04,  4.1839539D-03, -2.0193186D-02, -3.3875977D-03,
     X  -1.1682498D-03, -1.1402269D-03, -9.3456470D-05, -6.7233036D-05,
     X   9.0913315D-05, -8.9100615D-03,  9.6640034D-04,  3.1221683D-03,
     X   2.9887290D-04,  1.7246785D-04,  1.2439358D-04,  3.4472204D-06,
     X   9.3388499D-06, -4.3408535D-05,  8.6605648D-04,  9.5173471D-05,
     X   2.8956494D-02, -1.1006073D-04,  8.5819791D-04,  6.8065558D-05,
     X   4.3911027D-06,  1.0242931D-06, -5.6875754D-07,  3.4378192D-07,
     X   3.6377996D-06,  2.7727839D-05,  1.4606367D-05,  7.5029567D-03,
     X  -1.3647579D-03, -3.6079074D-03,  1.8793794D-05, -8.0081960D-05,
     X  -4.3339447D-06,  1.0184916D-06, -3.2200186D-07, -3.5035172D-08,
     X   5.7426186D-08, -4.2225760D-07, -2.4458611D-06,  1.3027750D-05,
     X  -5.9070403D-04,  1.2255553D-03, -3.9641933D-02, -5.2077168D-04,
     X  -5.1775053D-04, -2.9058748D-05, -1.2966185D-06, -3.5665825D-07/
C
      DATA (D_COEF(J,3,2,2),J=63,121)/  -1.3091110D-07, -6.5103016D-09,
     X   2.5653061D-09, -9.9290011D-09,  4.8959117D-09, -1.7042808D-07,
     X  -3.2562427D-06, -1.7409850D-05, -2.3560335D-04, -7.3923427D-03,
     X  -1.3669918D-02,  3.1584176D-03, -2.2343944D-04,  4.7205868D-05,
     X   1.7824652D-06, -1.4999564D-07,  5.2739462D-08,  2.4929673D-09,
     X   2.0320298D-09,  9.8075579D-11,  3.4677735D-10,  2.6475448D-09,
     X  -4.0253570D-09,  1.5420940D-07,  2.2899930D-06, -3.1548749D-06,
     X   4.9986959D-04, -1.6303059D-03,  4.7787282D-02,  7.5989869D-04,
     X   4.7140345D-04,  2.0902663D-05,  3.6608517D-07,  1.0216546D-07,
     X   9.8913820D-09, -4.0784565D-10,  1.4826661D-10, -8.8848049D-11,
     X  -8.1879968D-12, -3.6333526D-12, -6.4182010D-11,  1.2711921D-09,
     X  -5.4677351D-09,  3.0914733D-08,  6.5757825D-07,  1.0252274D-05,
     X   1.6523714D-04,  5.4006717D-03,  1.3714500D-02, -2.3500424D-03,
     X   1.1183319D-04, -2.4558721D-05, -8.2570420D-07,  2.1948788D-08,
     X  -3.9761688D-09, -1.5140502D-10, -4.5536779D-11, -1.5083076D-11,
     X  -3.0037083D-12/
C
      DATA (D_COEF(J,1,3,2),J=1,62)/     1.4307321D-03,  3.5272097D-03,
     X   3.0872513D-04, -5.1906120D-04, -1.0186446D-03,  6.0894669D-04,
     X  -5.5010597D-03, -4.9367560D-04, -1.1682202D-03, -1.8083704D-04,
     X   3.7320184D-05, -2.7429638D-03,  1.1318786D-03,  1.5136553D-03,
     X   4.8077831D-05,  1.8108404D-04,  1.8194809D-05, -2.4354427D-05,
     X   1.8670383D-04, -9.5097369D-04,  3.1365556D-03,  4.9284348D-04,
     X   1.5209539D-04,  5.0610442D-06, -1.8343050D-06,  1.0028661D-06,
     X   4.4849640D-06,  8.5152674D-07,  1.4112806D-04,  1.1922472D-03,
     X   4.3373787D-03, -1.0531248D-03,  1.0838185D-04, -1.9264255D-05,
     X   5.8555353D-07, -1.2249292D-06, -3.4008130D-08, -2.5442702D-07,
     X  -3.4476925D-06,  8.8962984D-06, -1.1212817D-04,  2.2372546D-03,
     X  -4.5761824D-03, -8.7873881D-04, -1.1943300D-04, -2.9812040D-05,
     X  -2.8167975D-06, -5.8779313D-08, -4.7936320D-08,  3.5952081D-09,
     X  -4.9661959D-09,  1.5549259D-07, -1.2520660D-06, -3.7145398D-06,
     X  -1.7953923D-04, -8.9610412D-04, -9.3533945D-03,  4.1569697D-04,
     X  -1.3909633D-04,  1.2500532D-05,  5.6214540D-07,  6.4988919D-08/
C
      DATA (D_COEF(J,1,3,2),J=63,121)/  -5.2826058D-09,  7.3071481D-09,
     X   1.8748805D-11, -3.5012055D-09,  4.5941883D-09, -1.0270858D-08,
     X   2.2358811D-07, -4.6705803D-06, -1.7386947D-05, -1.5664276D-03,
     X  -2.8720514D-03,  8.4701876D-04, -3.0069890D-05,  9.4528029D-06,
     X   3.9091297D-07,  5.0257693D-08, -2.7563644D-09, -1.4662718D-10,
     X   3.9134486D-10,  2.4745667D-11,  1.0935962D-12,  8.5810806D-10,
     X  -6.4108188D-09, -3.4000114D-08,  9.6900480D-08,  8.4611339D-07,
     X   2.9098302D-05, -3.1617772D-04,  3.6262885D-03,  2.6228265D-04,
     X   7.3266190D-06, -2.6758683D-06, -1.7256187D-08, -4.5435586D-08,
     X   4.7202124D-09,  6.6453598D-11,  1.4147205D-11, -6.1088624D-12,
     X  -1.9801134D-12,  4.9893674D-12, -5.9941106D-12,  2.2998922D-10,
     X  -1.8449246D-10, -9.5204814D-09,  1.8177565D-07, -7.5359318D-07,
     X   8.0027130D-06,  2.5111373D-05,  2.5200987D-03, -1.7716337D-04,
     X   3.5332691D-05, -1.9963288D-07,  4.1122342D-07, -1.3957429D-08,
     X  -2.6800696D-10,  1.4040768D-10,  2.9330808D-12, -7.9103415D-12,
     X  -1.0550335D-12/
C
      DATA (D_COEF(J,2,3,2),J=1,62)/     2.3375503D-04, -1.8790183D-03,
     X   9.0098955D-04,  3.6265738D-03,  4.7535418D-04,  8.6629801D-04,
     X   1.6825089D-03, -3.3457870D-04, -3.3868522D-05, -1.6023768D-04,
     X   1.9623938D-05,  1.0521483D-03, -3.0009149D-03, -5.2019082D-05,
     X  -3.6754958D-04, -2.6616210D-04, -9.0210360D-06, -3.3451019D-06,
     X  -2.2859982D-04, -1.0523683D-03, -3.0210772D-03,  1.1388260D-03,
     X   7.3113406D-05,  7.8359561D-05,  2.4952662D-05,  3.8375990D-06,
     X   2.3348419D-06, -2.4999321D-06,  1.0659562D-04, -2.8874885D-05,
     X   4.8573866D-03,  5.5505465D-04,  1.2364620D-04,  1.3384925D-05,
     X   2.7166238D-06, -1.5566290D-06, -1.1119749D-07, -4.1603579D-07,
     X   8.5901350D-07,  3.9880262D-07,  1.6501760D-04,  1.2114455D-03,
     X   1.3374418D-03, -5.0571656D-04,  1.0831197D-05, -1.0898693D-05,
     X  -5.3715275D-06,  2.3280334D-07, -3.3156750D-08,  3.7678886D-09,
     X  -1.5604510D-08, -1.4446061D-07, -6.3335251D-07, -3.5789234D-06,
     X  -1.0911237D-04,  2.2323667D-04, -5.6481840D-03, -5.4957512D-05,
     X  -4.7246161D-05, -9.7574104D-06, -3.2011575D-08, -2.5896732D-07/
C
      DATA (D_COEF(J,2,3,2),J=63,121)/  -1.0420483D-08,  9.6341906D-09,
     X  -3.3333545D-10, -1.2458429D-09,  2.0146240D-08,  1.0411423D-08,
     X  -7.3066553D-07, -1.0791219D-05, -3.7342812D-05, -1.3490182D-03,
     X  -1.6419033D-03,  3.1201934D-04, -2.2806578D-05,  8.1432120D-06,
     X   7.5654774D-07,  5.0054548D-08,  7.6042916D-09,  5.8859721D-10,
     X   2.1596725D-11,  7.0415632D-12,  5.4963206D-11,  1.6643306D-09,
     X  -1.6496582D-09,  3.2760214D-08,  2.4771346D-07, -1.4188590D-06,
     X   8.7445601D-05, -3.5315389D-04,  5.5803393D-03,  4.2292058D-05,
     X   2.5544975D-05,  1.7163877D-06, -6.0467665D-07,  2.7812043D-09,
     X   4.0161560D-09,  5.0731822D-10, -6.6316438D-11,  1.5974755D-11,
     X  -8.6615872D-13, -3.3710847D-12, -3.9438786D-11, -2.2205547D-10,
     X   1.5200598D-09, -2.8358023D-08,  3.7113899D-07, -6.5571507D-07,
     X  -4.4473995D-05,  1.1514715D-03,  2.0070093D-03, -7.0915636D-04,
     X   8.1293813D-05, -5.8032474D-06,  6.5359095D-07,  1.1555348D-08,
     X  -1.9661720D-09,  5.6373057D-10,  4.2143996D-11, -8.3193352D-12,
     X  -1.5004342D-12/
C
      DATA (D_COEF(J,3,3,2),J=1,62)/     1.1794120D-03,  4.4730621D-03,
     X  -5.2240925D-03, -2.1425612D-03, -3.0358405D-04, -5.8691155D-03,
     X  -1.7204988D-03,  1.2857678D-03, -4.6005429D-04, -2.7471574D-04,
     X   8.4443072D-04, -1.4280715D-03,  1.4278692D-02,  9.6846041D-04,
     X   9.8859698D-04,  3.8748559D-04,  4.9238519D-05,  6.5742713D-05,
     X   8.5146729D-05,  6.6472128D-03,  4.8199845D-04, -2.6041907D-03,
     X  -7.6075052D-06, -1.5502766D-04, -3.9959438D-05, -3.3234813D-06,
     X  -8.8209096D-06,  2.6640873D-05, -7.2793424D-04,  8.0709837D-04,
     X  -2.5331502D-02, -3.6163983D-04, -7.1139586D-04, -5.6616864D-05,
     X  -1.0290882D-06, -3.1504017D-06,  1.4570133D-07, -1.5745318D-07,
     X  -5.6408640D-06, -3.2816625D-05, -1.3286215D-04, -7.3025938D-03,
     X  -3.3224283D-03,  3.1658302D-03, -1.0217376D-04,  8.6571712D-05,
     X   2.8706205D-06, -6.0552309D-07,  4.6564073D-07,  3.6209128D-08,
     X  -5.2237543D-08,  4.8711662D-07,  3.0340399D-06, -7.6066546D-06,
     X   5.6872508D-04, -1.4987174D-03,  3.6821300D-02,  7.6046231D-04,
     X   5.2280768D-04,  2.9962301D-05,  3.2211506D-07,  3.4630447D-07/
C
      DATA (D_COEF(J,3,3,2),J=63,121)/   6.9648950D-08,  3.3480282D-09,
     X  -2.1518136D-09,  8.2023887D-09, -3.1149859D-09,  8.7851075D-08,
     X   2.4864414D-06,  1.7235329D-05,  2.1021923D-04,  6.7248225D-03,
     X   1.2214744D-02, -2.8925528D-03,  1.6862022D-04, -4.6708761D-05,
     X  -1.4673226D-06,  1.4845761D-07, -4.7123373D-08, -2.4565342D-10,
     X  -1.8204649D-09, -7.8482100D-11, -2.4188992D-10, -2.3639233D-09,
     X   5.9889750D-09, -1.3439658D-07, -2.0232026D-06,  6.3200212D-07,
     X  -4.4378420D-04,  1.3266835D-03, -4.2449379D-02, -7.7262156D-04,
     X  -4.1644193D-04, -1.6451432D-05,  1.1288669D-09, -8.3881165D-08,
     X  -7.7902750D-09,  3.3077772D-10, -2.2248858D-10,  8.0661817D-11,
     X   6.1884014D-12,  4.4091854D-12,  5.8662834D-11, -1.1234450D-09,
     X   4.4232778D-09, -2.0536567D-08, -3.2381773D-07, -8.3721445D-06,
     X  -1.2259999D-04, -4.5573363D-03, -1.0656259D-02,  2.0311518D-03,
     X  -6.7980838D-05,  2.1690111D-05,  5.9572149D-07, -2.8005340D-08,
     X   2.9083130D-09,  4.7714435D-11,  5.4007924D-11,  1.3308505D-11,
     X   3.1885939D-12/
c
c     Put the coefficients in COMMON so that they can be modified
c     by an initialization routine.  This allows the same program to
c     use different models for the magnetic coordinates.  The default
c     is defined by the coefficients in the DATA statements above.
c
      COMMON /sfc$$coeffs_com/D_COEF
C
C
C
C         The following IF-THEN-ELSE block checks the input arguments to
C         ensure they are within allowable ranges.  If any argument is
C         outside the acceptable range, the error flag, I_ERROR will be
C         set to some non-zero value.
C
      IF ((R_HEIGHT_IN .LT. 0.) .OR. (R_HEIGHT_IN .GT. 2000.)) THEN
C
C             The height in kilometers is outside the allowable range.
C
        I_ERROR = -2
C
      ELSEIF ((I_FLAG .LT. 1) .OR. (I_FLAG .GT. 2)) THEN
C
C             The conversion flag is neither 1 nor 2.
C
        I_ERROR = -4
C
      ELSEIF (ABS(R_LAT_IN) .GT. 90.) THEN
C
C             The latitude is outside the allowable range.
C
        I_ERROR = -8
C
      ELSEIF ((R_LON_IN .LT. 0.) .OR. (R_LON_IN .GT. 360.)) THEN
C
C             The longitude is outside the allowable range.
C
        I_ERROR = -16
C
      ELSE
C
C         All input arguments are within allowable ranges.
C
        I_ERROR = 0
C
      ENDIF
C
C
C         Quit if any error condition exists.
C
      IF (I_ERROR .NE. 0) RETURN
C
      IF (R_HEIGHT_IN .NE. R_HEIGHT_OLD(I_FLAG)) THEN
C
      ALT_VAR    = R_HEIGHT_IN/1200.0
      ALT_VAR_SQ = ALT_VAR * ALT_VAR
C
      DO 88 I = 1, 3

      DO 84 J = 1, 121
C
      D_CINT(J,I,I_FLAG) = D_COEF(J,I,1,I_FLAG) + ALT_VAR *
     x D_COEF(J,I,2,I_FLAG) + ALT_VAR_SQ *  D_COEF(J,I,3,I_FLAG)
C
 84   CONTINUE
C
 88   CONTINUE
C
        R_HEIGHT_OLD(I_FLAG) = R_HEIGHT_IN
C
      ENDIF
C
C     ZERO SUM FOR SPHERICAL HARMONIC EXPANSION COMPUTATION
C
      D_X = 0.0D0
      D_Y = 0.0D0
      D_Z = 0.0D0
C
C
CCCCCCCC ** PREPARE FOR COMPUTING SPHERICAL HARMONIC EXPANSION
C
C
      D_LON_INPUT   = DBLE(R_LON_IN) * DEGRAD
C
      IF (I_FLAG .EQ. 1) THEN
C
        D_COLAT_INPUT = (90.0D0 - DBLE(R_LAT_IN)) * DEGRAD
C
      ELSE
C
        CALL CG_ALT_DIP(R_HEIGHT_IN,R_LAT_IN,2,R_LAT_ADJ,IER64)
C
        D_COLAT_INPUT = (90.0D0 - DBLE(R_LAT_ADJ)) * DEGRAD
C
      ENDIF
C
CCCCCCCCC
C
C       Generate the spherical harmonics at the coordinate point
C
        CALL RYLM(D_COLAT_INPUT,D_LON_INPUT,I_ORDER,D_YLMVAL) 
C
C       Calculate the Cartesian coordinates of the unit vector
C       in the magnetic dipole coordinate system.
C
        DO 20 L= 0, I_ORDER
          DO 10 M = -L, L
C
            K = L * (L + 1) + M + 1
C
C           Add the contribution of each spherical harmonic to the
C           appropriate Cartesian component of the unit vector in the
C           magnetic dipole coordinate system.
C
C
            D_X = D_X + D_CINT(K,1,I_FLAG) * D_YLMVAL(K)
            D_Y = D_Y + D_CINT(K,2,I_FLAG) * D_YLMVAL(K)
            D_Z = D_Z + D_CINT(K,3,I_FLAG) * D_YLMVAL(K)
C
 10       CONTINUE
 20     CONTINUE
C
C           Compute the magnitude of the Cartesian unit vector of the
C           magnetic dipole coordinate system.
C
        D_R = SQRT(D_X * D_X + D_Y * D_Y + D_Z * D_Z)
C
C           If the magnitude of the unit vector differs significantly
C           from 1, set the error flag and continue processing.
C
        IF ((D_R .LT. 0.9D0) .OR. (D_R .GT. 1.1D0)) THEN
          I_ERROR = -32
        ENDIF
C
C           Adjust the components so they do represent the components of
C           a unit vector.  If D_R is equal to 1.0D0, this step will not
C           change the values of D_X, D_Y, or D_Z.
C
        D_Z = D_Z / D_R
        D_X = D_X / D_R
        D_Y = D_Y / D_R
C
C           Obtain output co_latitude and longitude from the unit
C           vector using standard formulas
C
        IF (D_Z .GT. 1.0D0) THEN
          D_COLAT_TEMP = 0.0D0
        ELSE IF (D_Z .LT.- 1.0D0) THEN
          D_COLAT_TEMP = PI
        ELSE
          D_COLAT_TEMP = DACOS(D_Z)
        ENDIF
C
        IF ((ABS(D_X) .LT. 1.0D-8) .AND. (ABS(D_Y) .LT. 1.0D-8)) THEN
          D_LON_TEMP = 0.0D0
        ELSE
          D_LON_TEMP = DATAN2(D_Y,D_X)
        ENDIF
C
CCCCCCCCC *** PREPARE LATITUDE DATA FOR OUTPUT
C
       D_LON_OUTPUT   = D_LON_TEMP
C
       IF (I_FLAG .EQ. 1) THEN
C
         R_LAT_TMP = 90.0 - D_COLAT_TEMP/DEGRAD
         CALL CG_ALT_DIP(R_HEIGHT_IN,R_LAT_TMP,1,R_LAT_ADJ,IER64)
C
         D_COLAT_OUTPUT = (90.0D0 - DBLE(R_LAT_ADJ)) * DEGRAD
C
       ELSE
C
         D_COLAT_OUTPUT = D_COLAT_TEMP
C
       ENDIF
C
CCCCCCCCC
C
C       Convert colatitude into latitude and convert both coordinates
C       from DOUBLE PRECISION radians to REAL degrees.
C
        R_LAT_OUT = SNGL(90.0D0 - D_COLAT_OUTPUT / DEGRAD)
        R_LON_OUT = SNGL(D_LON_OUTPUT / DEGRAD)
C
      IF (R_LON_OUT .LT. 0.0) R_LON_OUT = R_LON_OUT + 360.
C
      IF (IER64 .EQ. 1) I_ERROR = -64
      RETURN
      END

