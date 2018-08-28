#include <stdio.h>

void cnv_mdhms_sec_idl(argc,argv)
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

  *seconds = cnv_mdhms_sec(yr,mo,dy,hr,mn,sc);

}
