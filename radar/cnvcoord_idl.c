#include <stdio.h>

void cnvcoord_idl(argc,argv)
     int argc;
     void *argv[];
{
  short *mgflag;
  short *order;
  short *err;
  float *inpos;
  float *in_lat;
  float *in_long;
  float *height;
  float *outpos;
  float *out_lat;
  float *out_long;
  float *out_r;

  inpos  = (float *) argv[0];
  order  = (short *) argv[1];
  outpos = (float *) argv[2];
  mgflag = (short *) argv[3];
  err    = (short *) argv[4];

  in_lat  = inpos;
  in_long = (inpos+1);
  height  = (inpos+2);

  out_lat  = outpos;
  out_long = (outpos+1);
  out_r    = (outpos+2);
/*  outpos = inpos; */
 
 cnv$coord_(in_lat,in_long,height,order,out_lat,out_long,out_r,mgflag,err);
   
}
