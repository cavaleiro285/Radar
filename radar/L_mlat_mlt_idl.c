#include <stdio.h>

void L_mlat_mlt_idl(argc,argv)

    int argc;
    void *argv[];
{
    short   *iyr, *iday, *ierr;
    long  *isod;
    float dum,*r, *theta, *phi, *l, *mlat, *mlon, *mlt, *outpos;

    dum = 0.0;
    iyr   = (short *) argv[0];
    iday  = (short *) argv[1];
    isod  = (long *) argv[2];
    r     = (float *) argv[3];
    theta = (float *) argv[4];
    phi   = (float *) argv[5];
    outpos = (float *) argv[6]; 
    ierr  = (short *) argv[7];

    l = &dum;
    mlon = outpos;
    mlat = (outpos+1);
    mlt = (outpos+2);
    

#ifdef __ALPHA
     l_mlat_mlt_(iyr,iday,isod,r,theta,phi,l,mlat,mlon,mlt,ierr);
#endif

#ifdef SUNOS
     _l_mlat_mlt_(iyr,iday,isod,r,theta,phi,l,mlat,mlon,mlt,ierr);
#endif

#ifdef AIX
   l_mlat_mlt(iyr,iday,isod,r,theta,phi,l,mlat,mlon,mlt,ierr);
#endif

#ifdef HPUX
  l_mlat_mlt(iyr,iday,isod,r,theta,phi,l,mlat,mlon,mlt,ierr);
#endif




}
