c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license
c
c ig1 contains the column numbers of group 1 (others are assumed to be in
c group2).  ig2 is an integer working array
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

c on output dist(1)=prop permutations with at least nlr(1) stats<=crit(1)
c dist(2)=prop permutations with any stats<=crit(1)
c dist(3)=ave # stats <=crit(1) (dist(4:6) similarly for upper tail).
c values of crit determined from ordered stats and values of nlr
      subroutine ptn(d,n,ng,ng1,z,pzil,pzir,pzml,pzmr,nperm,ix,
     $     nclust,mclust,mct1,ig1,ig2,irnk,istrt,wght,
     $     nlr,crit,dist,iflag,ipc,zt)
      real d(n,ng),z(n),ztml,ztmr,tsum,wght(nclust),crit(2)
      real dist(6),zt(n)
      integer pzml(n),pzmr(n)
      integer nclust, mclust(nclust), mct1(nclust), ig1(ng),ig2(ng)
      integer irnk,istrt,pzil(n),pzir(n),iflag(nclust),ipc,nct1
      integer n,ng,ng1,i,j,k,l,nperm,isamp,ix(3),inl,inr,nlr(2)
      if (ipc.eq.1) then
         nct1=0
         do 2 i=1,nclust
            if (mct1(i).gt.0) then
               nct1=nct1+1
c next line only relevant if nperm=0
               ig2(nct1)=i
            endif
 2       continue 
      endif
      do 3 i=1,6
         dist(i)=0
 3    continue 
      if (nperm.le.0) then
         isamp=0
         nperm=1
         dist(1)=1
         dist(2)=1
         dist(3)=nlr(1)
         dist(4)=1
         dist(5)=1
         dist(6)=nlr(2)
         do 16 i=1,n
            pzil(i)=1
            pzir(i)=1
            pzml(i)=1
            pzmr(i)=1
 16      continue 
      else
         isamp=1
         dist(1)=0
         dist(2)=0
         dist(3)=0
         dist(4)=0
         dist(5)=0
         dist(6)=0
         do 17 i=1,n
            pzil(i)=0
            pzir(i)=0
            pzml(i)=0
            pzmr(i)=0
 17      continue 
      endif
 20   if (isamp.gt.0) then
         if (isamp.gt.nperm) go to 80
         isamp=isamp+1
         if (ipc.ne.1) then
            call rpg(ng1,ng,nclust,mclust,mct1,ig1,ig2,ix)
         else
c choose nct1 clusters for group 1      
            call rnsub(nct1,nclust,ig2,ig1,ix)
            call sortg(ig2,ig1,1,nct1)
c ng1 can change in this call
            call pgc(ng1,ng,nclust,mclust,nct1,ig1,ig2)
         endif
      else
         if (ipc.ne.1) then
            k=1
            l=1
            do 410 j=1,nclust
               call upc(ig1(l),mct1(j),k+mclust(j)-1,k,iflag(j))
               if (iflag(j).lt.1) then
                  do 411 i=1,j-1
                     iflag(i)=0
 411              continue 
                  go to 27
               endif
               k=k+mclust(j)
               l=l+mct1(j)
 410        continue 
         else
            call upc(ig2,nct1,nclust,1,iflag)
c ng1 can change in this call
            call pgc(ng1,ng,nclust,mclust,nct1,ig1,ig2)
            if (iflag(1).lt.1) go to 27
         endif
c no more permutations
         go to 80
 27      nperm=nperm+1
      endif
      inl=0
      inr=0
      ztmr=z(1)-1
      do 40 i=1,n
         zt(i)=tsum(ng1,d(i,1),n,ig1,istrt,nclust,mclust,mct1,wght)
         if (zt(i).ge.z(i)) pzir(i)=pzir(i)+1
         if (zt(i).le.z(i)) pzil(i)=pzil(i)+1
         if (zt(i).le.crit(1)) inl=inl+1
         if (zt(i).ge.crit(2)) inr=inr+1
         ztmr=max(zt(i),ztmr)
         if (ztmr.ge.z(i)) pzmr(i)=pzmr(i)+1
 40   continue
      ztml=z(n)+1
      do 43 i=n,1,-1
         ztml=min(zt(i),ztml)
         if (ztml.le.z(i)) pzml(i)=pzml(i)+1
 43   continue 
      if (inl.ge.nlr(1)) dist(1)=dist(1)+1
      if (inl.ge.1) dist(2)=dist(2)+1
      dist(3)=dist(3)+inl
      if (inr.ge.nlr(2)) dist(4)=dist(4)+1
      if (inr.ge.1) dist(5)=dist(5)+1
      dist(6)=dist(6)+inr
      go to 20
 80   continue
c enforce monotonicity
      do 45 i=2,n
         pzml(i)=max(pzml(i-1),pzml(i))
         pzmr(n-i+1)=max(pzmr(n-i+1),pzmr(n-i+2))
 45   continue 
      do 82 i=1,6
         dist(i)=dist(i)/nperm
 82   continue
      return
      end
