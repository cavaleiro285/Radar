	subroutine l_mlat_mlt(syr,sday,isod,r,
     1		theta,phi,l,mlat,mlon,mlt,serr)
c
c	Input parameters:
c		iyr =	Year to use for the calculation
c		iday = 	Day of year -- used for MLT computation
c    		ihr =	Hour of day -- used for MLT computation
c		imin =	Minute of hour -- used for MLT computation
c		isec =	Second of minute -- used for MLT computation
c    NOTE       isod =  Seconds of day -- replaces ihr,imin,isec
c		r   =	Radial distance of point at which to calculate
c		theta =	Latitude of point (in radians)
c		phi =	Longitude of point (in radians)
c
c	Output parameters:
c		l =	L-shell parameter
c		mlat =	Magnetic Latitude
c               mlon =  Magnetic Longitude
c		mlt =	Magnetic Local Time
c		ierr =	Error flag
c
	implicit none
c
        real pi,degrad,radeg 
	real r,theta,phi,l,mlat,mlon,mlt,g(66),h(66),x,y,z
	real alpha,beta,x0,y0,z0,xd,yd,zd,xn,yn,zn,re
	real ca,sa,cb,sb,thetam,cthm,cthm2,phim,phimd
	real gst,slong,srasn,sdec,sx,sy,sz,sxn,syn,szn,phis
	real sgst,cgst,sxgeo,sygeo,szgeo
	integer*2 syr,sday,serr
	integer*4 iyr,iday,isod,ierr
c
        pi=4.*atan(1.0)
        degrad=pi/180.
        radeg=180./pi
c patch

        iyr=int(syr)
        iday=int(sday)
c       ihr=int(shr)
c       imin=int(smin)
c       isec=int(ssec)
 	ierr=int(serr)
        theta=theta*degrad
        phi=phi*degrad
c
	call get_igrf_coeffs(iyr,g,h,ierr)
c
c	Convert input point from (r,theta,phi) to GEO (x,y,z)
c
	x=r*cos(theta)*cos(phi)
	y=r*cos(theta)*sin(phi)
	z=r*sin(theta)
c
c	Compute terms for rotation matrix and shift to
c	eccentric dipole coordinates. These terms depend
c	on the IGRF coefficients, and thus on which 
c	specific IGRF model is selected.
c
	call get_terms(g,h,alpha,beta,x0,y0,z0)
c
c	Displace origin to eccentric dipole coordinates (xd,yd,zd)
c
	re=6378.160
c	write(3,*) 'x0,y0,z0 = ',x0*re,y0*re,z0*re
c
	xd = x - x0
	yd = y - y0
	zd = z - z0
c
c	Rotate to system aligned with diple in azimuth and tilt
c
c	[xn]=[Cos(alpha)Cos(beta)  Cos(alpha)Sin(beta)  -Sin(alpha)] [xd]
c	[yn]=[    -Sin(beta)             Cos(beta)          0      ] [yd]
c	[zn]=[Sin(alpha)Cos(beta)  Sin(alpha)Sin(beta)   Cos(alpha)] [zd]
c
	ca = cos(alpha)
	sa = sin(alpha)
	cb = cos(beta)
	sb = sin(beta)
c
	xn = ca*cb*xd + ca*sb*yd - sa*zd
	yn = -sb*xd + cb*yd
	zn = sa*cb*xd + sa*sb*yd + ca*zd
c
c	Compute position of the sun -- needed to compute MLT
c
c        print *, iyr,iday,isod
	call sun(iyr,iday,isod,gst,slong,srasn,sdec)
c
 	sx = cos(srasn)*cos(sdec)
 	sy = sin(srasn)*cos(sdec)
	sz = sin(sdec)
c
c	Transform to GEO coordinates (rotation about z-axis by
c	angle corresponding to mean Grenwich sideral time).
c
	sgst=sin(gst)
	cgst=cos(gst)

	sxgeo = cgst*sx + sgst*sy
	sygeo = -sgst*sx + cgst*sy
	szgeo = sz
c
c	Rotate sun position into system aligned with dipole tilt in 
c	azimuth and tilt -- i.e. multiply by the same matrix the we
c	used above to transform (xd,yd,zd) into (xn,yn,zn).
c
        sxn = ca*cb*sxgeo + ca*sb*sygeo - sa*szgeo
        syn = -sb*sxgeo + cb*sygeo
        szn = sa*cb*sxgeo + sa*sb*sygeo + ca*szgeo
c
	phis = atan2(syn,sxn)
c
c	Compute magnetic latitude
c
	thetam = atan(zn/sqrt(xn*xn+yn*yn))
	mlat = thetam*radeg
c
c	Compute L-shell parameter
c
	cthm = cos(thetam)
	cthm2 = cthm*cthm
	l = r/cthm2
c
c	Compute magnetic local time (note addition of pi to make
c	MLT=0 to correpond to midnight and MLT=12 to noon)
c
	phim = atan2(yn,xn)
	phimd = (phim - phis + pi)*radeg
        mlon = phim*radeg
c
c	write(*,*) 'Diagnostics:'
c	write(*,*)
c	write(*,*) 'xd,yd,zd = ',xd,yd,zd
c	write(*,*) 'xn,yn,zn = ',xn,yn,zn
c	write(*,*) 'sxgeo,sygeo,szgeo = ',sxgeo,sygeo,szgeo
c	write(*,*) 'sxn,syn,szn = ',sxn,syn,szn
c	write(*,*) 'phim, phis = ',phim*radeg,phis*radeg
c	write(*,*) 'phimd = ',phimd
c	write(*,*)
c
c	Force phimd to be in (0-360) range, then divide by 15 to
c	convert to hours in (0-24) range
c
	if (phimd.lt.0.0) phimd = phimd + 360.
	if (phimd.gt.360.) phimd = phimd - 360.
	mlt = phimd/15.

c        print *, mlon,mlat,mlt

	return
	end


C
	subroutine get_igrf_coeffs(iyr,g,h,ierr)
c
C	SET UP TO ACCEPT DATES BETWEEN 1965 AND 2000; COEFFICIENTS
C	THROUGH 1985 ARE FROM DGRF MODELS COEFFICIENTS FOR 1990 FROM IGRF
C	[EOS TRANS. AGU APRIL 21, 1992, P. 182]. INTERPOLATION IS USED
c	FOR YEARS BETWEEN DGRF MODELS, AND EXTRAPOLATION FOR YEARS
C	BETWEEN 1990 AND 2000.
c
c
c		Mauricio Peredo
c		Hughes STX at NASA/GSFC
c		June 29, 1995
c
C  Additional IGRF References:
C  (J. GEOMAG. GEOELECTR.(1982), V.34, P.313-315,
C  GEOMAGN. AND AERONOMY (1986), V.26, P.523-525).
C
C History
C
C DGRF90 & IGRF 95 added.  R. Baldwin HSTX   9/95
C
C------INPUT PARAMETERS:
C  iyr - YEAR NUMBER (FROM 1965 UP TO 2000)
C----- OUTPUT PARAMETERS:
C  g,h - coefficients for the igrf model interpolated (or extrapolated)
c	 to the epoch of iyr
C
C
	IMPLICIT NONE
C
      REAL F1,F2,DT,G(66),H(66),G65(66),H65(66),G70(66),H70(66),
     *G75(66),H75(66),G80(66),H80(66),G85(66),H85(66),G90(66),
     *H90(66),G95(66),H95(66),G00(66),H00(66),DG00(45),DH00(45)

	INTEGER iyr,N,ierr
c
      DATA G65/0.,-30334.,-2119.,-1662.,2997.,1594.,1297.,-2038.,1292.,
     *856.,957.,804.,479.,-390.,252.,-219.,358.,254.,-31.,-157.,-62.,
     *45.,61.,8.,-228.,4.,1.,-111.,75.,-57.,4.,13.,-26.,-6.,13.,1.,13.,
     *5.,-4.,-14.,0.,8.,-1.,11.,4.,8.,10.,2.,-13.,10.,-1.,-1.,5.,1.,-2.,
     *-2.,-3.,2.,-5.,-2.,4.,4.,0.,2.,2.,0./
      DATA H65/0.,0.,5776.,0.,-2016.,114.,0.,-404.,240.,-165.,0.,148.,
     *-269.,13.,-269.,0.,19.,128.,-126.,-97.,81.,0.,-11.,100.,68.,-32.,
     *-8.,-7.,0.,-61.,-27.,-2.,6.,26.,-23.,-12.,0.,7.,-12.,9.,-16.,4.,
     *24.,-3.,-17.,0.,-22.,15.,7.,-4.,-5.,10.,10.,-4.,1.,0.,2.,1.,2.,
     *6.,-4.,0.,-2.,3.,0.,-6./
c
      DATA G70/0.,-30220.,-2068.,-1781.,3000.,1611.,1287.,-2091.,1278.,
     *838.,952.,800.,461.,-395.,234.,-216.,359.,262.,-42.,-160.,-56.,
     *43.,64.,15.,-212.,2.,3.,-112.,72.,-57.,1.,14.,-22.,-2.,13.,-2.,
     *14.,6.,-2.,-13.,-3.,5.,0.,11.,3.,8.,10.,2.,-12.,10.,-1.,0.,3.,
     *1.,-1.,-3.,-3.,2.,-5.,-1.,6.,4.,1.,0.,3.,-1./
      DATA H70/0.,0.,5737.,0.,-2047.,25.,0.,-366.,251.,-196.,0.,167.,
     *-266.,26.,-279.,0.,26.,139.,-139.,-91.,83.,0.,-12.,100.,72.,-37.,
     *-6.,1.,0.,-70.,-27.,-4.,8.,23.,-23.,-11.,0.,7.,-15.,6.,-17.,6.,
     *21.,-6.,-16.,0.,-21.,16.,6.,-4.,-5.,10.,11.,-2.,1.,0.,1.,1.,3.,
     *4.,-4.,0.,-1.,3.,1.,-4./
c
      DATA G75/0.,-30100.,-2013.,-1902.,3010.,1632.,1276.,-2144.,1260.,
     *830.,946.,791.,438.,-405.,216.,-218.,356.,264.,-59.,-159.,-49.,
     *45.,66.,28.,-198.,1.,6.,-111.,71.,-56.,1.,16.,-14.,0.,12.,-5.,
     *14.,6.,-1.,-12.,-8.,4.,0.,10.,1.,7.,10.,2.,-12.,10.,-1.,-1.,4.,
     *1.,-2.,-3.,-3.,2.,-5.,-2.,5.,4.,1.,0.,3.,-1./
      DATA H75/0.,0.,5675.,0.,-2067.,-68.,0.,-333.,262.,-223.,0.,191.,
     *-265.,39.,-288.,0.,31.,148.,-152.,-83.,88.,0.,-13.,99.,75.,-41.,
     *-4.,11.,0.,-77.,-26.,-5.,10.,22.,-23.,-12.,0.,6.,-16.,4.,-19.,6.,
     *18.,-10.,-17.,0.,-21.,16.,7.,-4.,-5.,10.,11.,-3.,1.,0.,1.,1.,3.,
     *4.,-4.,-1.,-1.,3.,1.,-5./
c
      DATA G80/0.,-29992.,-1956.,-1997.,3027.,1663.,1281.,-2180.,1251.,
     *833.,938.,782.,398.,-419.,199.,-218.,357.,261.,-74.,-162.,-48.,
     *48.,66.,42.,-192.,4.,14.,-108.,72.,-59.,2.,21.,-12.,1.,11.,-2.,
     *18.,6.,0.,-11.,-7.,4.,3.,6.,-1.,5.,10.,1.,-12.,9.,-3.,-1.,7.,2.,
     *-5.,-4.,-4.,2.,-5.,-2.,5.,3.,1.,2.,3.,0./
      DATA H80/0.,0.,5604.,0.,-2129.,-200.,0.,-336.,271.,-252.,0.,212.,
     *-257.,53.,-297.,0.,46.,150.,-151.,-78.,92.,0.,-15.,93.,71.,-43.,
     *-2.,17.,0.,-82.,-27.,-5.,16.,18.,-23.,-10.,0.,7.,-18.,4.,-22.,9.,
     *16.,-13.,-15.,0.,-21.,16.,9.,-5.,-6.,9.,10.,-6.,2.,0.,1.,0.,3.,
     *6.,-4.,0.,-1.,4.,0.,-6./
c
      DATA G85/0.,-29873.,-1905.,-2072.,3044.,1687.,1296.,-2208.,1247.,
     *829.,936.,780.,361.,-424.,170.,-214.,355.,253.,-93.,-164.,-46.,
     *53.,65.,51.,-185.,4.,16.,-102.,74.,-62.,3.,24.,-6.,4.,10.,0.,21.,
     *6.,0.,-11.,-9.,4.,4.,4.,-4.,5.,10.,1.,-12.,9.,-3.,-1.,7.,1.,-5.,
     *-4.,-4.,3.,-5.,-2.,5.,3.,1.,2.,3.,0./
      DATA H85/0.,0.,5500.,0.,-2197.,-306.,0.,-310.,284.,-297.,0.,232.,
     *-249.,69.,-297.,0.,47.,150.,-154.,-75.,95.,0.,-16.,88.,69.,-48.,
     *-1.,21.,0.,-83.,-27.,-2.,20.,17.,-23.,-7.,0.,8.,-19.,5.,-23.,11.,
     *14.,-15.,-11.,0.,-21.,15.,9.,-6.,-6.,9.,9.,-7.,2.,0.,1.,0.,3.,
     *6.,-4.,0.,-1.,4.,0.,-6./
c  new
      DATA G90/0., -29775.,  -1848.,  -2131.,   3059.,   1686.,   1314.,
     *     -2239.,   1248.,    802.,    939.,    780.,    325.,   -423.,
     *       141.,   -214.,    353.,    245.,   -109.,   -165.,    -36.,
     *        61.,     65.,     59.,   -178.,      3.,     18.,    -96.,
     *        77.,    -64.,      2.,     26.,     -1.,      5.,      9.,
     *         0.,     23.,      5.,     -1.,    -10.,    -12.,      3.,
     *         4.,      2.,     -6.,      4.,      9.,      1.,    -12.,
     *         9.,     -4.,     -2.,      7.,      1.,     -6.,     -3.,
     *        -4.,      2.,     -5.,     -2.,      4.,      3.,      1.,
     *         3.,      3.,      0./
  
      DATA H90/0.,      0.,   5406.,      0.,  -2279.,   -373.,      0.,
     *      -284.,    293.,   -352.,      0.,    247.,   -240.,     84.,
     *      -299.,      0.,     46.,    154.,   -153.,    -69.,     97.,
     *         0.,    -16.,     82.,     69.,    -52.,      1.,     24.,
     *         0.,    -80.,    -26.,      0.,     21.,     17.,    -23.,
     *        -4.,      0.,     10.,    -19.,      6.,    -22.,     12.,
     *        12.,    -16.,    -10.,      0.,    -20.,     15.,     11.,
     *        -7.,     -7.,      9.,      8.,     -7.,      2.,      0.,
     *         2.,      1.,      3.,      6.,     -4.,      0.,     -2.,
     *         3.,     -1.,     -6./

      DATA G95/0., -29682.,  -1789.,  -2197.,   3074.,   1685.,   1329.,
     *     -2268.,   1249.,    769.,    941.,    782.,    291.,   -421.,
     *       116.,   -210.,    352.,    237.,   -122.,   -167.,    -26.,
     *        66.,     64.,     65.,   -172.,      2.,     17.,    -94.,
     *        78.,    -67.,      1.,     29.,      4.,      8.,     10.,
     *        -2.,     24.,      4.,     -1.,     -9.,    -14.,      4.,
     *         5.,      0.,     -7.,      4.,      9.,      1.,    -12.,
     *         9.,     -4.,     -2.,      7.,      0.,     -6.,     -3.,
     *        -4.,      2.,     -5.,     -2.,      4.,      3.,      1.,
     *         3.,      3.,      0./
  
      DATA H95/0.,      0.,   5318.,      0.,  -2356.,   -425.,      0.,
     *      -263.,    302.,   -406.,      0.,    262.,   -232.,     98.,
     *      -301.,      0.,     44.,    157.,   -152.,    -64.,     99.,
     *         0.,    -16.,     77.,     67.,    -57.,      4.,     28.,
     *         0.,    -77.,    -25.,      3.,     22.,     16.,    -23.,
     *        -3.,      0.,     12.,    -20.,      7.,    -21.,     12.,
     *        10.,    -17.,    -10.,      0.,    -19.,     15.,     11.,
     *        -7.,     -7.,      9.,      7.,     -8.,      1.,      0.,
     *         2.,      1.,      3.,      6.,     -4.,      0.,     -2.,
     *         3.,     -1.,     -6./

      DATA G00/0., -29615.,  -1728.,  -2267.,   3072.,   1672.,   1341.,
     *     -2290.,   1253.,    715.,    935.,    787.,    251.,   -405.,
     *       110.,   -217.,    351.,    222.,   -131.,   -169.,    -12.,
     *        72.,     68.,     74.,   -161.,     -5.,     17.,    -91.,
     *        79.,    -74.,      0.,     33.,      9.,      7.,      8.,
     *        -2.,     25.,      6.,     -9.,     -8.,    -17.,      9.,
     *         7.,     -8.,     -7.,      5.,      9.,      3.,     -8.,
     *         6.,     -9.,     -2.,      9.,     -4.,     -8.,     -2.,
     *        -6.,      2.,     -3.,      0.,      4.,      1.,      2.,
     *         4.,      0.,     -1./
  
      DATA H00/0.,      0.,   5186.,      0.,  -2478.,   -458.,      0.,
     *      -227.,    296.,   -492.,      0.,    272.,   -232.,    119.,
     *      -304.,      0.,     44.,    172.,   -134.,    -40.,    107.,
     *         0.,    -17.,     64.,     65.,    -61.,      1.,     44.,
     *         0.,    -65.,    -24.,      6.,     24.,     15.,    -25.,
     *        -6.,      0.,     12.,    -22.,      8.,    -21.,     15.,
     *         9.,    -16.,     -3.,      0.,    -20.,     13.,     12.,
     *        -6.,     -8.,      9.,      4.,     -8.,      5.,      0.,
     *         1.,      0.,      4.,      5.,     -6.,     -1.,     -3.,
     *         0.,     -2.,     -8./

      DATA DG00/0.0,  14.6,    10.7,   -12.4,     1.1,    -1.1,     0.7,
     *         -5.4,   0.9,    -7.7,    -1.3,     1.6,    -7.3,     2.9,
     *         -3.2,   0.0,    -0.7,    -2.1,    -2.8,    -0.8,     2.5,
     *          1.0,  -0.4,     0.9,     2.0,    -0.6,    -0.3,     1.2,
     *         -0.4,  -0.4,    -0.3,     1.1,     1.1,    -0.2,     0.6,
     *         -0.9,  -0.3,     0.2,    -0.3,     0.4,    -1.0,     0.3,
     *         -0.5,  -0.7,    -0.4/
  
      DATA DH00/0.0,   0.0,   -22.5,     0.0,   -20.6,    -9.6,     0.0,
     *          6.0,  -0.1,   -14.2,     0.0,     2.1,     1.3,     5.0,
     *          0.3,   0.0,    -0.1,     0.6,     1.7,     1.9,     0.1,
     *          0.0,  -0.2,    -1.4,     0.0,    -0.8,     0.0,     0.9,
     *          0.0,   1.1,     0.0,     0.3,    -0.1,    -0.6,    -0.7,
     *          0.2,   0.0,     0.1,     0.0,     0.0,     0.3,     0.6,
     *         -0.4,   0.3,     0.7/


c      DATA DG95/0.0,  17.6,    13.0,   -13.2,     3.7,    -0.8,     1.5,
c     *         -6.4,  -0.2,    -8.1,     0.8,     0.9,    -6.9,     0.5,
c     *         -4.6,   0.8,     0.1,    -1.5,    -2.0,    -0.1,     2.3,
c     *          0.5,  -0.4,     0.6,     1.9,    -0.2,    -0.2,     0.0,
c     *         -0.2,  -0.8,    -0.6,     0.6,     1.2,     0.1,     0.2,
c     *         -0.6,   0.3,    -0.2,     0.1,     0.4,    -1.1,     0.3,
c     *          0.2,  -0.9,    -0.3/
c  
c      DATA DH95/0.0,   0.0,   -18.3,     0.0,   -15.0,    -8.8,     0.0,
c     *          4.1,   2.2,   -12.1,     0.0,     1.8,     1.2,     2.7,
c     *         -1.0,   0.0,     0.2,     1.2,     0.3,     1.8,     0.9,
c     *          0.0,   0.3,    -1.6,    -0.2,    -0.9,     1.0,     2.2,
c     *          0.0,   0.8,     0.2,     0.6,    -0.4,     0.0,    -0.3,
c     *          0.0,   0.0,     0.4,    -0.2,     0.2,     0.7,     0.0,
c     *         -1.2,  -0.7,    -0.6/
c  old
c     DATA G90/0.,-29775.,-1851.,-2136.,3058.,1693.,1315.,-2240.,1246.,
c    *807.,939.,782.,324.,-423.,142.,-211.,353.,244.,-111.,-166.,-37.,
c    *61.,64.,60.,-178.,2.,17.,-96.,77.,-64.,4.,28.,1.,6.,10.,0.,22.,
c    *5.,-1.,-11.,-12.,4.,4.,3.,-6.,4.,10.,1.,-12.,9.,-4.,-1.,7.,2.,-6.,
c    *-4.,-4.,2.,-5.,-2.,4.,3.,1.,2.,3.,0./
c     DATA H90/0.,0.,5411.,0.,-2278.,-380.,0.,-287.,293.,-348.,0.,248.,
c    *-240.,87.,-299.,0.,47.,153.,-154.,-69.,98.,0.,-16.,83.,68.,-52.,
c    *2.,27.,0.,-81.,-27.,1.,20.,16.,-23.,-5.,0.,10.,-20.,7.,-22.,12.,
c    *6.,-4.,0.,-1.,4.,0.,-6./
c     DATA DG90 /0.,18.0,10.6,-12.9,2.4,0.0,3.3,-6.7,0.1,-5.9,0.5,0.6,
c    *-7.0,0.5,-5.5,0.6,-0.1,-1.6,-3.1,-0.1,2.3,1.3,-0.2,1.8,1.3,-0.2,
c    *0.1,1.2,0.6,-0.5,-0.3,0.6,1.6,0.2,0.2,0.3,0.2,-0.7,-0.2,0.1,-1.1,
c    *0.0,-0.1,-0.5,-0.6/
c     DATA DH90 /0.,0.,-16.1,0.,-15.8,-13.8,0.,4.4,1.6,-10.6,0.,2.6,1.8,
c    *3.1,-1.4,0.,-0.1,0.5,0.4,1.7,0.4,0.,0.2,-1.3,0.0,-0.9,0.5,1.2,
c    *0.,0.6,0.2,0.8,-0.5,-0.2,0.0,0.0,0.,0.5,-0.2,0.3,0.3,0.4,-0.5,
c    *-0.3,0.6/
c
C
	if ( (iyr.lt.1965).or.(iyr.gt.2005) ) then
	   ierr=1
	   write(*,999) iyr
	   goto 300
	endif
c
      IF (IYR.LT.1970) GOTO 50		!INTERPOLATE BETWEEN 1965 - 1970
      IF (IYR.LT.1975) GOTO 60		!INTERPOLATE BETWEEN 1970 - 1975
      IF (IYR.LT.1980) GOTO 70		!INTERPOLATE BETWEEN 1975 - 1980
      IF (IYR.LT.1985) GOTO 80		!INTERPOLATE BETWEEN 1980 - 1985
      IF (IYR.LT.1990) GOTO 90		!INTERPOLATE BETWEEN 1985 - 1990
      IF (IYR.LT.1995) GOTO 100		!INTERPOLATE BETWEEN 1990 - 1995
      IF (IYR.LT.2000) GOTO 110		!INTERPOLATE BETWEEN 1995 - 2000
C
C	EXTRAPOLATE BETWEEN 2000 - 2005
C
      DT=FLOAT(IYR)-2000.
      DO 40 N=1,66
         G(N)=G00(N)
         H(N)=H00(N)
         IF (N.GT.45) GOTO 40
         G(N)=G(N)+DG00(N)*DT
         H(N)=H(N)+DH00(N)*DT
40    CONTINUE
      GOTO 300
C
C
C	INTERPOLATE BETWEEEN 1965 - 1970
C
50    F2=(IYR-1965)/5.
      F1=1.-F2
      DO 55 N=1,66
         G(N)=G65(N)*F1+G70(N)*F2
55       H(N)=H65(N)*F1+H70(N)*F2
      GOTO 300
C
C
C	INTERPOLATE BETWEEN 1970 - 1975
C
60    F2=(IYR-1970)/5.
      F1=1.-F2
      DO 65 N=1,66
         G(N)=G70(N)*F1+G75(N)*F2
65       H(N)=H70(N)*F1+H75(N)*F2
      GOTO 300
C
C
C	INTERPOLATE BETWEEN 1975 - 1980
C
70    F2=(IYR-1975)/5.
      F1=1.-F2
      DO 75 N=1,66
         G(N)=G75(N)*F1+G80(N)*F2
75       H(N)=H75(N)*F1+H80(N)*F2
      GOTO 300
C
C
C	INTERPOLATE BETWEEN 1980 - 1985
C
80    F2=(IYR-1980)/5.
      F1=1.-F2
      DO 85 N=1,66
         G(N)=G80(N)*F1+G85(N)*F2
85       H(N)=H80(N)*F1+H85(N)*F2
      GOTO 300
C
C
C	INTERPOLATE BETWEEN 1985 - 1990
C
90    F2=(IYR-1985)/5.
      F1=1.-F2
      DO 95 N=1,66
         G(N)=G85(N)*F1+G90(N)*F2
95       H(N)=H85(N)*F1+H90(N)*F2
      GOTO 300
C
c
c      Interpolate between 1990 - 1995
C
100   F2=(IYR-1990)/5.
      F1=1.-F2
      DO 105 N=1,66
       G(N)=G90(N)*F1 + G95(N)*F2
105    H(N)=H90(N)*F1 + H95(N)*F2
      GOTO 300       
C
c
c      Interpolate between 1995 - 2000
C
110   F2=(IYR-1995)/5.
      F1=1.-F2
      DO 115 N=1,66
       G(N)=G95(N)*F1 + G00(N)*F2
115    H(N)=H95(N)*F1 + H00(N)*F2
      GOTO 300       
C
300	continue
c
C
999	format(//1x,
     *	 '*** ERROR -- Input year = ',i5,/
     *   ' is out of valid range 1965-2005 ***'//)
c
	return
	end


C
c
	subroutine get_terms(g,h,alpha,beta,x0,y0,z0)
c
	implicit none
c
	real g(66),h(66),alpha,beta,x0,y0,z0
	real g10,g11,h11,g20,g21,g22,h21,h22
	real b02,b0,sq3,l0,l1,l2,e
c
	g10 = g(2)
	g11 = g(3)
	h11 = h(3)
	g20 = g(4)
	g21 = g(5)
	g22 = g(6)
	h21 = h(5)
	h22 = h(6)
c
	b02 = g10*g10 + g11*g11 + h11*h11
	b0 = sqrt(b02)
	alpha = acos(-g10/b0)
	beta = atan(h11/g11)
c
	sq3 = sqrt(3.)
	l0 = 2.*g10*g20 + sq3*(g11*g21 + h11*h21)
	l1 = -g11*g20 + sq3*(g10*g21 + g11*g22 + h11*h22)
	l2 = -h11*g20 + sq3*(g10*h21 - h11*g22 + g11*h22)
	e = (l0*g10 + l1*g11 + l2*h11)/(4.*b02)
c
	z0 = (l0 - g10*e)/(3.*b02)
	x0 = (l1 - g11*e)/(3.*b02)
	y0 = (l2 - h11*e)/(3.*b02)
c
	return
	end

c--------------------------------------------------------------------
c
      SUBROUTINE SUN(IYR,IDAY,ISOD,GST,SLONG,SRASN,SDEC)
C
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  THIS SUBROUTINE HAS BEEN COMPILED FROM: RUSSELL C.T., COSM.ELECTRO-
C  DYN., 1971, V.2,PP.184-196.
C
C
C                   AUTHOR: Gilbert D. Mead
C
C
        IMPLICIT NONE

        REAL GST,SLONG,SRASN,SDEC,RAD,T,VL,G,OBLIQ,SOB,SLP,SIND,
     1       COSD,SC

        INTEGER IYR,IDAY,ISOD

      DOUBLE PRECISION DJ,FDAY
      DATA RAD/57.295779513/
      IF(IYR.LT.1901.OR.IYR.GT.2099) RETURN
c      FDAY=DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      FDAY=DFLOAT(ISOD)/86400.D0
      DJ=365*(IYR-1900)+(IYR-1901)/4+IDAY-0.5D0+FDAY
      T=DJ/36525.
      VL=DMOD(279.696678+0.9856473354*DJ,360.D0)
      GST=DMOD(279.690983+.9856473354*DJ+360.*FDAY+180.,360.D0)/RAD
      G=DMOD(358.475845+0.985600267*DJ,360.D0)/RAD
      SLONG=(VL+(1.91946-0.004789*T)*SIN(G)+0.020094*SIN(2.*G))/RAD
      IF(SLONG.GT.6.2831853) SLONG=SLONG-6.2831853
      IF (SLONG.LT.0.) SLONG=SLONG+6.2831853
      OBLIQ=(23.45229-0.0130125*T)/RAD
      SOB=SIN(OBLIQ)
      SLP=SLONG-9.924E-5
C
C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
C   THE ORBITAL MOTION OF THE EARTH
C
      SIND=SOB*SIN(SLP)
      COSD=SQRT(1.-SIND**2)
      SC=SIND/COSD
      SDEC=ATAN(SC)
      SRASN=3.141592654-ATAN2(COS(OBLIQ)/SOB*SC,-COS(SLP)/COSD)
      RETURN
      END
c-----------------------------------------------------------------------
c  




