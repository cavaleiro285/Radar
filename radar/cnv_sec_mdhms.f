       subroutine cnv_sec_mdhms(yr,mo,dy,hr,min,sc,t0)
       integer*2 yr,mo,dy,hr,min,sc
       integer*4 modays(12) 
       data modAys/31,28,31,30,31,30,31,31,30,31,30,31/
       integer*4 t0
        
       if(mod(yr,4).eq.0) modays(2)=29 

       iday = t0/86400
       nday = iday*86400
       nsec = t0 - nday
       ihr = nsec/3600
       nhr = ihr*3600
       nsec = nsec - nhr
       imin = nsec/60
       nmin = imin*60
       isec = nsec - nmin
       if(iday.eq.0) iday=1
       
       isum =0
       last =0
       do i=1,12
        if(isum.lt.iday) then
         isum=isum+modays(i)
         last=modays(i)
         imo=i
        end if
       end do

       iday=iday-(isum-last) 

       mo=imo
       dy=iday
       hr=ihr
       min=imin
       sc=isec
       
       return
       end
