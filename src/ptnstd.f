c By R Gray, 
c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license
c
c ptnstd:
c ig1 contains the column numbers of group 1 (others are assumed to be in
c group2).  
c If nclust>1, then columns of d must be sorted on cluster, with the first
c    mclust(1) columns in cluster 1, the next mclust(2) in cluster 2, etc.
c    mct1(j) is the number of columns in group 1 in cluster j
c if istrt=1 then tests are stratified on the groups defined by cluster
c    in this case standardization is done within stratum, and the statistic
c    is computed as the sum over strata of the group 1 means of the
c    standardized values within each stratum
c if ipc=1 then clusters are permuted among groups, rather than columns 
c    within clusters 
c Note: single precision stats
      subroutine ptnstd(d,n,ng,ng1,z,nclust,mclust,mct1,ig1,irnk,istrt,
     $     wght,ipc)
      real d(n,ng),z(n),tsum,wght(nclust)
      integer nclust, mclust(nclust), mct1(nclust), ig1(ng1)
      integer irnk,istrt,ipc
      integer n,ng,ng1,i,j,k
      if (istrt.eq.1) then
c standardize within cluster
         k=1
         do 320 j=1,nclust
c adjust weights so a weight of 1 gives within strata difference of means
            if (mct1(j).gt.0.and.mct1(j).lt.mclust(j)) then
               wght(j)=wght(j)*mclust(j)/(mct1(j)*(mclust(j)-mct1(j)))
            else
               wght(j)=0
            endif
            if (irnk.eq.1) then
c means only (in a sense variances already stratified since values are ranks)
               call stdm(n,d(1,k),mclust(j))
            else if (irnk.ne.2) then
c means and se
               call stdmv(n,d(1,k),mclust(j))
            endif
            k=k+mclust(j)
 320     continue 
      else
         if (irnk.eq.1) then
c means only
            call stdm(n,d,ng)
         else if (irnk.ne.2) then
c means and se
            call stdmv(n,d,ng)
         endif
      endif
      do 10 i=1,n
         z(i)=tsum(ng1,d(i,1),n,ig1,istrt,nclust,mclust,mct1,wght)
 10   continue 
      return
      end
