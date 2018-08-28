#include <stdio.h>

void cnv_sec_mdhms_idl(argc,argv)
     int argc;
     void *argv[];
{
  short *yr, *mo, *dy, *hr, *mn, *sc;
  long *seconds;

  yr = (short *) argv[0];
  mo = (short *) argv[1];
  dy = (short *) argv[2];
  hr = (short *) argv[3];
  mn = (short *) argv[4];
  sc = (short *) argv[5];

  seconds = (long *) argv[6];

/*  cnv_sec_mdhms_(yr,mo,dy,hr,mn,sc,seconds);
 */
#ifdef __ALPHA
    cnv_sec_mdhms_(yr,mo,dy,hr,mn,sc,seconds);
#endif

#ifdef SUNOS
     _cnv_sec_mdhms_(yr,mo,dy,hr,mn,sc,seconds);
#endif

#ifdef AIX
    cnv_sec_mdhms(yr,mo,dy,hr,mn,sc,seconds);
#endif

#ifdef HPUX
    cnv_sec_mdhms(yr,mo,dy,hr,mn,sc,seconds);
#endif

}
