#include <stdio.h>
#include <math.h>

float mlt_idl(argc,argv)
     int argc;
     void *argv[];
{
  short *year;
  long  *t;
  float *mlong,mslong,mt;

#ifdef SUNOS
  extern FLOATFUNCTIONTYPE mlt_();
#else
  extern float mlt_(),_mlt_(),mlt();
#endif


  year  = (short *) argv[0];
  t     = (long * ) argv[1];
  mlong = (float *) argv[2];


#ifdef __ALPHA 
     mt = mlt_(year,t,mlong,&mslong);
#endif

#ifdef SUNOS
     ASSIGNFLOAT(mt, mlt_(year,t,mlong,&mslong));
#endif

#ifdef AIX
    mt = mlt(year,t,mlong,&mslong);
#endif

#ifdef HPUX 
    mt = mlt(year,t,mlong,&mslong);
#endif

  return mt;
}

