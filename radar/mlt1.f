	real*4 function MLT1(T0,solar_dec,mlong1,mslong)
c
c     $Revision: 1.1 $
c
c     $Log:	mlt1.f,v $
c Revision 1.1  94/10/14  11:28:49  11:28:49  baker (Kile Baker S1G)
c Initial revision
c 
c
	implicit none
	real*8 T0,told,t2
	real*4 solar_dec,mlong1,mslong,mslong1,mslong2,mslat1,mslat2
	real*4 slong,slong1,slong2,mrad,mslat,sol_dec_old,height
	integer*2 ier,order,mflag
        integer*2 itemp
        data itemp/0/
	save

c
c	if the time hasn't changed by more than 10 minutes
c	interpolate the sun's magnetic position.
c

	if ((abs(solar_dec - sol_dec_old) .gt. .01)
     &		.or. (sol_dec_old .eq. 0.0)) told = 86400  !force new calc
	if (abs(mslong2 - mslong1).gt. 10.0) told = 86400
	if ((t0 .ge. told) .and. (t0 .lt. t2)) then
	  mslong = mslong1 + (t0-told)*(mslong2-mslong1)/(t2-told)  !interpolate
	else
	  told = t0
	  sol_dec_old = solar_dec
	  slong1 = (12.*3600. - t0)*15./3600.
	  t2 = t0 + 600				!t0 plus 10 minutes
	  slong2 = (12.*3600. - t2)*15./3600.
          height = 450.0
          order = 4
          mflag = 1
	  call cnv$coord(solar_dec,slong1,height,order
     &         ,mslat1,mslong1,mrad,mflag,ier)
	  call cnv$coord(solar_dec,slong2,height,order
     &         ,mslat2,mslong2,mrad,mflag,ier)
	  mslong = mslong1
	end if
	MLT1 = (mlong1 - mslong)/15.0 + 12.0
	if (mlt1.ge.24.0) mlt1 = mlt1 - 24.0
	if (mlt1.lt.0.) mlt1 = mlt1 + 24.0
	return
	end
