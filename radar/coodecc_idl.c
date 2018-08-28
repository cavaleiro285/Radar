#include <stdio.h>

void coodecc_idl(argc,argv)

int argc;
void *argv[];

{
   short *iyr, *irot, *ierr;
   float *lat, *lon, *r;

   lat = (float *) argv[0];
   lon = (float *) argv[1];
   r   = (float *) argv[2];
   iyr = (short *) argv[3];
   irot= (short *) argv[4];
   ierr= (short *) argv[5];

#ifdef __ALPHA
     coodecc_(lat,lon,r,iyr,irot,ierr);
#endif

#ifdef SUNOS
     coodecc_(lat,lon,r,iyr,irot,ierr);
#endif

#ifdef AIX
    coodecc(lat,lon,r,iyr,irot,ierr);
#endif

#ifdef HPUX
    coodecc(lat,lon,r,iyr,irot,ierr);
#endif

}
