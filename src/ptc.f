c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

c counts complete separation
c if istrt=1 only makes comparisons within clusters
c if nclust>1 then only permutes within clusters
c ig1 gives the column numbers of subjects in group 1
c ig2 gives the column numbers of subjects in group 2
c Note: single precision and integers only
      subroutine ptc(d,n,ng,ng1,ics,nperm,dtcs,ix,nclust,mclust,mct1,
     $     ig1,ig2,istrt,iflag,ipc,igc)
      real d(n,ng)
      integer dtcs(4),istrt,ipc,igc(nclust),nct1,iflag(nclust)
      integer nclust, mclust(nclust),mct1(nclust),ig1(ng1),ig2(ng)
      integer ix(3),nperm,it,ics(n)
c      integer mg,
c      parameter(mg=200)
      integer n,ng,ng1,i,j,k,l,ng2,isamp,icmps
c initialize
      if (ipc.eq.1) then
         nct1=0
         do 2 i=1,nclust
            if (mct1(i).gt.0) then
               nct1=nct1+1
c next line only relevant if nperm=0
               igc(nct1)=i
            endif
 2       continue 
      endif
      ng2=ng-ng1
      dtcs(1)=0
      do 10 i=1,n
         call tstatc(d(i,1),ig1,ng1,ig2,ng2,n,icmps,nclust,mclust,mct1,
     $        istrt)
         ics(i)=icmps
         dtcs(1)=dtcs(1)+icmps
 10   continue 
      if (nperm.le.0) then
         isamp=0
         nperm=1
         dtcs(2)=1
         dtcs(3)=dtcs(1)
         dtcs(4)=min(dtcs(1),1)
      else
         isamp=1
         dtcs(2)=0
         dtcs(3)=0
         dtcs(4)=0
      endif
 20   if (isamp.gt.0) then
         if (isamp.gt.nperm) go to 80
         isamp=isamp+1
         if (ipc.ne.1) then
            call rpg(ng1,ng,nclust,mclust,mct1,ig1,ig2,ix)
            go to 810
         else
c choose nct1 clusters for group 1      
            call rnsub(nct1,nclust,igc,ig1,ix)
            call sortg(igc,ig1,1,nct1)
c ng1 can change in this call
            call pgc(ng1,ng,nclust,mclust,nct1,ig1,igc)
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
            call upc(igc,nct1,nclust,1,iflag)
c ng1 can change in this call
            call pgc(ng1,ng,nclust,mclust,nct1,ig1,igc)
            if (iflag(1).lt.1) go to 27
         endif
c no more permutations
         go to 80
 27      nperm=nperm+1
      endif
c skip this if random permutations and ipc=0
      ng2=ng-ng1
      k=1
      do 30 j=1,ng
         if (k.le.ng1) then
            if (j.eq.ig1(k)) then
               k=k+1
               go to 30
            endif
         endif
         ig2(j-k+1)=j
 30   continue 
 810  it=0
      do 40 i=1,n
         call tstatc(d(i,1),ig1,ng1,ig2,ng2,n,icmps,nclust,mclust,mct1,
     $        istrt)
         it=it+icmps
 40   continue 
      if (it.ge.dtcs(1)) dtcs(2)=dtcs(2)+1
      dtcs(3)=dtcs(3)+it
      if (it.gt.0) dtcs(4)=dtcs(4)+1
      go to 20
 80   return
      end
