c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

      function tsum(ng,d,n,ig1,istrt,nclust,mclust,mct1,wght)
      real d(n,*),tsum,wght(nclust)
      integer n,ng,ig1(ng),i,j,l,istrt,nclust,mclust(nclust)
      integer mct1(nclust)
      double precision dd,dd1
      dd=0
      if (istrt.ne.1) then
         do 10 i=1,ng
            dd=dd+d(1,ig1(i))
 10      continue 
      else
         l=0
         do 12 j=1,nclust
            if (mct1(j).gt.0.and.mct1(j).lt.mclust(j)) then
               dd1=0
               do 13 i=1,mct1(j)
                  dd1=dd1+d(1,ig1(l+i))
 13            continue 
c note that data was standardized and the weights adjusted in ptnstd so that 
c wght(j)*dd1 is the input weight to ptn times the difference of the means
               dd=dd+wght(j)*dd1
            endif
            l=l+mct1(j)
 12      continue
      endif
      tsum=dd
      return
      end
