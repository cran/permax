c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

      subroutine stdm(n,d,ng)
      integer n,ng,i,j
      real d(n,ng)
      double precision dmm
      do 10 i=1,n
         dmm=0
         do 12 j=1,ng
            dmm=dmm+d(i,j)
 12      continue 
         dmm=dmm/ng
         do 14 j=1,ng
            d(i,j)=(d(i,j)-dmm)
 14      continue 
 10   continue 
      return
      end

      subroutine stdmv(n,d,ng)
      integer n,ng,i,j
      real d(n,ng)
      double precision dmm,dss
      do 10 i=1,n
         dmm=0
         do 12 j=1,ng
            dmm=dmm+d(i,j)
 12      continue 
         dmm=dmm/ng
         dss=0
         do 13 j=1,ng
c do not use d(i,j)=d(i,j)-dmm here, because it converts to single precision
            dss=dss+(d(i,j)-dmm)**2
 13      continue 
         dss=sqrt(dss/(ng-1))
         if (dss.gt.0) then
            do 14 j=1,ng
               d(i,j)=(d(i,j)-dmm)/dss
 14         continue 
         endif
 10   continue 
      return
      end

c for correlations -- don't divide by ng-1
      subroutine stdms(n,d,ng)
      integer n,ng,i,j
      real d(n,ng)
      double precision dmm,dss
      do 10 i=1,n
         dmm=0
         do 12 j=1,ng
            dmm=dmm+d(i,j)
 12      continue 
         dmm=dmm/ng
         dss=0
         do 13 j=1,ng
            dss=dss+(d(i,j)-dmm)**2
 13      continue 
         dss=sqrt(dss)
         if (dss.gt.0) then
            do 14 j=1,ng
               d(i,j)=(d(i,j)-dmm)/dss
 14         continue 
         endif
 10   continue 
      return
      end
