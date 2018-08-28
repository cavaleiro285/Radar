#include <math.h>
 
cnv_sec_mdhms(iyr,mo,dy,hr,min,sec,seconds)
long *seconds;
short *iyr,*mo,*dy,*hr,*min,*sec;
{
double fmod(),year,myear,dsec,spm,mph,hpd;
double dmin,dhr,ddy;
short days,hours,tdys;
short modys[12];
short leap,i;

modys[0]=modys[2]=modys[4]=modys[6]=modys[7]=modys[9]=modys[11]=31;
modys[3]=modys[5]=modys[8]=modys[10]=30;
modys[1] = 28;

myear = 4.;
year = (double)*iyr;
if( fmod(year,myear) == 0. && *mo > 2) modys[1]=29;
dsec=(double)*seconds;
spm=60.;
mph=60.;
hpd=24.;
*sec=(short)fmod(dsec,spm);
dmin=(double)(*seconds-*sec)/60.;
*min=(short)fmod(dmin,mph);
dhr=(dmin-(double)*min)/60.;
*hr=(short)fmod(dhr,hpd);
ddy=(dhr-(double)*hr)/24.;
*mo = (short)(ddy/365.*12.)+1;
tdys=0;
for( i=0; i<*mo-1; i++) tdys += modys[i];
if( ddy-tdys < 0 ){
  *mo--;
  tdys=0;
  for( i=0; i<*mo-1; i++) tdys += modys[i];
}
*dy = ddy-tdys+1;
if( *dy > modys[*mo-1] ){
  *dy = 1;
  ++*mo;
}
}

