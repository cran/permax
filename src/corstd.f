c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

      subroutine corstd(d,n,ng,phen,z,nclust,mclust,istrt,wght)
      real dip,d(n,ng),phen(ng),z(n),wght(nclust)
      integer n,ng,i,j,k
      integer nclust,mclust(nclust),istrt
      if (istrt.eq.1) then
c standardize within stratum
         k=1
         do 320 j=1,nclust
c means and ss
            call stdms(n,d(1,k),mclust(j))
            call stdms(1,phen(k),mclust(j))
            k=k+mclust(j)
 320     continue 
      else
c means and ss
         call stdms(n,d,ng)
         call stdms(1,phen,ng)
      endif
c compute correlations
      do 311 i=1,n
         z(i)=dip(ng,d(i,1),n,phen,1,nclust,mclust,istrt,wght)
 311  continue
      return
      end
