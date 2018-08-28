	real*4 function mlt(iyr,t0,mlong,mslong)
c
c       $Revision: 1.2 $
c
c       This routine implements the calculation of Magnetic Local Time
c
c       It uses the Altitude Adjusted Corrected Geomag. Coordinates to
c       define the magnetic coordinates of the sun and the location of
c       the point where MLT is desired.
c
c       It also uses the equation of time from the Navy Ephemeris to
c       accurately calculate the positon of the sun.
c
c       The Inputs are: 
c           iyr   (integer*2)  the year (may be in the form nn or 19nn or 20nn)
c           t0    (integer*4)  the time in seconds since the beginning of year
c           mlong (real*4)     the magnetic longitude of the observation pt.
c
c       The outputs are:
c           the function returns a real*4 value giving the MLT time in hours
c           the variable "mslong" returns the magnetic longitude of the sun.
c
c       functions called:
c           solar_loc, eqn_of_time, cnv_sec_mdhms, mlt1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       $Log:	mlt.f,v $
c Revision 1.2  94/10/17  13:27:16  13:27:16  baker (Kile Baker S1G)
c added the year to the call to eqn_of_time
c 
c Revision 1.1  94/10/14  11:28:29  11:28:29  baker (Kile Baker S1G)
c Initial revision
c 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	integer*2 iyr,yr
	integer*4 t0,t1
	real*4 t2
	equivalence (t1,t2)
	real*4 mlong,mslong
	real*4 mean_long,dec
	integer*2 yrday,mo,day,hr,minut,sc
	real*8 ut,et,eqn_of_time,apparent_time
	real*4 mlt1
	save

	yr = iyr
	if (iyr.gt.1900) yr = iyr -1900
	if ((t0 .lt. 0.0) .or. (t0 .gt. 86400*366)) then !time error: its probably floating point
	  t1=t0				!move variable to local storage
	  t0=t2				!force conversion from real*4 to Int*4
	  type *,'MLT-E-UTFLOAT: argument 2 should be integer*4'
	  type *,'MLT-I-UTCONVRT: argument 2 reassigned value ',t0
	end if
	call solar_loc(yr,t0,mean_long,dec)
        
	et = eqn_of_time(mean_long, iyr)
        
	call cnv_sec_mdhms(yr,mo,day,hr,minut,sc,t0)
       
	ut = (hr*60.0 + minut)*60. + sc
	apparent_time = ut + et
	mlt = mlt1(apparent_time,dec,mlong,mslong)
        
	return
	end
