c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

c randomly permute u
c ix=seeds for Wichman Hill generator
      subroutine rnperm(ng,u,ix,nclust,mclust)
      integer ng,ix(3),i,ii,j,k,nclust,mclust(nclust),nx
      real whu,u(ng),tmp
      if (nclust.le.1) then
c at each step randomly select the ng-i+1st value from the set that has
c not been selected yet (in 1 : (ng-i+1))
         do 20 i=1,ng-1
            ii=whu(ix)*(ng-i+1)+1
            if (ii.lt.1.or.ii.gt.ng-i+1) write (6,*) 'error in rnperm'
            tmp=u(ng-i+1)
            u(ng-i+1)=u(ii)
            u(ii)=tmp
 20      continue 
      else
         k=1
         do 10 j=1,nclust
            nx=mclust(j)
            do 15 i=1,nx
               ii=whu(ix)*(nx-i+1)+k
               if (ii.lt.k.or.ii.gt.nx-i+k) write (6,*)'error in rnperm'
               tmp=u(nx-i+k)
               u(nx-i+k)=u(ii)
               u(ii)=tmp
 15         continue 
            k=k+mclust(j)
 10      continue 
      endif
      return
      end
