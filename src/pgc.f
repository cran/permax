c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

      subroutine pgc(m,ng,nclust,mclust,nct1,ig1,ig2)
      integer ng,nclust,mclust(nclust),nct1,ig1(ng),ig2(ng)
      integer i,j,k,l,m
      k=0
      l=1
      m=0
      do 10 i=1,nclust
         if (i.eq.ig2(l)) then
            do 11 j=1,mclust(i)
               m=m+1
               ig1(m)=k+j
 11         continue 
            if (l.eq.nct1) return
            l=l+1
         endif
         k=k+mclust(i)
 10   continue 
      return
      end
