      subroutine coodecc(lat,lon,r,iyr,irot,ierr)
      implicit none
      real lat,lon,r,rmag(3),rmag1(3),g(66),h(66),geo2ecc(3,3),
     *     ecc2geo(3,3),x0,y0,z0,dr,dlat,dlon,vm,vs,alpha,beta 
      integer iyr,irot,ierr

      data dr/0.017453293/

c Form rotation matrix for eccentric dipole
      
      call get_igrf_coeffs(iyr,g,h,ierr)                                      
 
c       Compute terms for rotation matrix and shift to
c       eccentric dipole coordinates. These terms depend
c       on the IGRF coefficients, and thus on which
c       specific IGRF model is selected.
c
      call get_terms(g,h,alpha,beta,x0,y0,z0)
c
c       Rotate to system aligned with diple in azimuth and tilt
c
c       [xn]=[Cos(alpha)Cos(beta)  Cos(alpha)Sin(beta)  -Sin(alpha)] [xd]
c       [yn]=[    -Sin(beta)             Cos(beta)          0      ] [yd]
c       [zn]=[Sin(alpha)Cos(beta)  Sin(alpha)Sin(beta)   Cos(alpha)] [zd]
c
      geo2ecc(1,1) = cos(alpha)*cos(beta)
      geo2ecc(1,2) = cos(alpha)*sin(beta)
      geo2ecc(1,3) = -sin(alpha)
      geo2ecc(2,1) = -sin(beta)
      geo2ecc(2,2) = cos(beta)
      geo2ecc(2,3) = 0.0
      geo2ecc(3,1) = sin(alpha)*cos(beta)
      geo2ecc(3,2) = sin(alpha)*sin(beta)
      geo2ecc(3,3) = cos(alpha)

C compute transpose 

      call tmat(ecc2geo,geo2ecc,3)

c Compute cartestian position from polar position
      lat = 90.0 - lat ! Convert latitude to colatitude
      if(lon.lt.0.0) lon=lon+360.0
      rmag(1)=r*sin(lat*dr)*cos(lon*dr)
      rmag(2)=r*sin(lat*dr)*sin(lon*dr)
      rmag(3)=r*cos(lat*dr)
 
c Perform the rotation based on the variable irot

      if(irot.eq.0) then

        rmag(1) = rmag(1) - x0
        rmag(2) = rmag(2) - y0
        rmag(3) = rmag(3) - z0

        call multmat(geo2ecc,rmag,rmag1)

      else   

        call multmat(ecc2geo,rmag,rmag1)

        rmag(1) = rmag(1) + x0
        rmag(2) = rmag(2) + y0
        rmag(3) = rmag(3) + z0

      endif 

        vm=sqrt(rmag1(1)**2+rmag1(2)**2+rmag1(3)**2)
        vs=rmag1(3)/vm
        dlat=asin(vs)
        dlon=atan2(rmag1(2),rmag1(1))
        lat=dlat/dr
        lon=dlon/dr
        r=vm

      return 
      end

      subroutine tmat(at,a,n)
      real*4 at(3,3),a(3,3)
      do i=1,n
       do j=1,n
        at(j,i)=a(i,j)
       end do
      end do
      return
      end
      SUBROUTINE MULTMAT (MatA, MatB, MatOut)

      REAL*4 MatA(3,3),                     ! 3x3 matrix
     +       MatB(3),                       ! 3x1 matrix
     +       MatOut(3)                      ! 3x1 result matrix
C*********************************************************************
C     Matrix multiplication                                          *
C                                                                    *
C        | Out1 |     | A11  A12  A13 |   | B1 |                     *
C        | Out2 |  =  | A21  A22  A23 | * | B2 |                     *
C        | Out3 |     | A31  A32  A33 |   | B3 |                     *
C                                                                    *
C*********************************************************************
      DO I = 1, 3
        MatOut(I) = MatA(i,1) * MatB(1) +
     1              MatA(i,2) * MatB(2) +
     2              MatA(i,3) * MatB(3)
      END DO
      RETURN
      END
