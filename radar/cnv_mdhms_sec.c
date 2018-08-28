/*	This routine converts the date (month, day, hour, minute, second)
	into seconds from the beginning of the year.  The year is required
	in order to make the adjustment for leap-years. */
/*
$Log: cnv_mdhms_sec.c,v $
 * Revision 1.4  1994/07/11  16:11:18  bristow
 * Changed argument declaratinion to short
 *
 * Revision 1.3  1994/07/08  11:08:08  baker
 * Changed the arguments from value to pointer.  This allows fortran
 * ro;utines to call this function.  This is different from the way
 * the routine is written for QNX.
 *
 * Revision 1.2  1993/06/15  20:24:38  kile
 * corrected the name of the revision variable
 *
 * Revision 1.1  1993/06/08  17:45:13  kile
 * Initial revision
 *
*/

char cnv_mdhms_sec_rev[] = {"$Revision: 1.4 $"};

/* static global definition of day of the year for each month (day 1) */
int jday[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
int mday[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

long cnv_mdhms_sec(yr, mo, day, hr, mn, sec)
	short *yr, *mo, *day, *hr, *mn, *sec;
{
long t;

/* first check for valid input values */

if ( (*mo < 1) || (*mo > 12) || (*hr < 0) || (*hr > 23) || 
        (*day <1) || (*mn < 0) ||
	(*mn > 59) || (*sec < 0) || (*sec >= 60) )
	        return -1;

if (*day > mday[*mo-1])
	if (*mo != 2)
		return -1;
	else if (*day != (mday[1] +1) || ((*yr % 4) != 0))
		return -1;

/*  ok, the input parameters are valid */

t = jday[*mo-1] + *day - 1;

if ( (*mo > 2) && ((*yr % 4) == 0)) t = t+1; /* inc day if leap year and mo>Feb */

return (((t*24L + *hr)*60L + *mn)*60L + *sec);
}

