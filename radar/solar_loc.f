	subroutine solar_loc(yr,t1,mean_long,declination)
	implicit none
c
c     $Revision: 1.2 $
c
c
c	This subroutine converts the year and time of the year (from
c	cnv$mdhms_sec) to the mean solar longitude of the sun and the
c	solar declination.
c
c*****************************************************************
c     $Log:	solar_loc.f,v $
c Revision 1.2  94/10/17  12:37:51  12:37:51  baker (Kile Baker S1G)
c corrected some of the constants
c 
c Revision 1.1  94/10/14  11:29:42  11:29:42  baker (Kile Baker S1G)
c Initial revision
c 
c
c*****************************************************************
	integer*2 yr,yr_step,delta_yr,i,iyr
	integer*4 t1
	integer*4 index
	real*4 mean_long,declination
	real*8 d,lambda,g,eps,L
	real*8 L0(10)/279.642,279.403,279.165,278.926,279.673,279.434,
     1                279.196,278.957,279.704,279.465/
	real*8 DL/0.985647/
	real*8 G0(10)/356.892984,356.637087,356.381191,356.125295,
     1                356.854999,356.599102,356.343206,356.087308,
     2                356.817011,356.561113/
	real*8 DG/0.98560028/
	real*8 EPS0(10)/23.440722,23.440592,23.440462,23.440332,
	1    23.440202,23.440072,23.439942,23.439811,23.439679,23.439548/
	real*8 DE/-0.00000036/
	save

	d = 0

	if (yr .lt. 1900) then
	   index = yr - 88
	else
	   index = yr - 1988
	endif
c	print *,'input values = ',yr,t1
c	print *, 'index = ',index
	if (index .le. 0) then
	   delta_yr = index - 1
	else if (index .gt. 10) then
	   delta_yr = index - 10
	else 
	   delta_yr = 0
	endif

	if (index .le. 0) index = 1
	if (index .gt. 10) index = 10

	yr_step = sign(1,delta_yr)
	delta_yr = abs(delta_yr)

c	print *,'delta_yr =',delta_yr
	do i = 1,delta_yr
	  if (yr_step .gt. 0) then
	    iyr = 98 + i - 1
	  else 
	    iyr = 89 - i
	  end if
	  if (mod(iyr,4).eq.0)then
	    d = d + 366*yr_step
	  else
	    d = d + 365*yr_step
	  end if
	end do
c	print *,'before using t1, d=',d
	d = d + t1/86400.
c	print *,'after using t1, d=',d
	L = L0(index) + DL*d
	g = G0(index) + DG*d
	do while (L.lt.0)
	  L = L + 360.
	end do
	do while (g.lt.0)
	  g = g + 360.
	end do

	L = mod(L,360.0d0)
	g = mod(g,360.0d0)

	lambda = L + 1.915*sind(g) + 0.020*sind(2*g)
	eps = EPS0(index) + DE*d

	declination = asind(sind(eps)*sind(lambda))
	mean_long = L
	return
	end

