c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

c determine whether this sample has complete separation
      subroutine tstatc(d,ig1,ng1,ig2,ng2,n,icmps,nclust,mclust,mct1,
     $     istrt)
      real d(n,ng1+ng2),di1,da1,di2,da2
      integer ig1(ng1),ng1,ig2(ng2),ng2,n,i,icmps,istrt,nclust,j,k,l
      integer mclust(nclust),mct1(nclust),idir
      if (istrt.ne.1) then
         di1=d(1,ig1(1))
         da1=di1
         di2=d(1,ig2(1))
         da2=di2
         do 10 i=2,ng1
            di1=min(di1,d(1,ig1(i)))
            da1=max(da1,d(1,ig1(i)))
 10      continue 
         do 12 i=2,ng2
            di2=min(di2,d(1,ig2(i)))
            da2=max(da2,d(1,ig2(i)))
 12      continue 
         if (da1.lt.di2.or.di1.gt.da2) then
            icmps=1
         else
            icmps=0
         endif
      else
c only check within clusters
         icmps=0
         idir=0
         l=0
         k=0
         do 112 j=1,nclust
            if (mct1(j).le.0.or.mct1(j).ge.mclust(j)) go to 113
            di1=d(1,ig1(l+1))
            da1=di1
            di2=d(1,ig2(k+1))
            da2=di2
            do 210 i=l+2,l+mct1(j)
               di1=min(di1,d(1,ig1(i)))
               da1=max(da1,d(1,ig1(i)))
 210        continue 
            do 212 i=k+2,k+mclust(j)-mct1(j)
               di2=min(di2,d(1,ig2(i)))
               da2=max(da2,d(1,ig2(i)))
 212        continue
c separation must be in the same direction in each cluster
            if (idir.eq.0) then
               if (di1.gt.da2) then
                  idir=1
               else if (da1.lt.di2) then
                  idir=2
               else
                  return
               endif
            else if (idir.eq.1) then
               if (di1.le.da2) return
            else 
               if (da1.ge.di2) return
            endif
 113        l=l+mct1(j)
            k=k+mclust(j)-mct1(j)
 112     continue 
         icmps=1
      endif
      return
      end
