c By R Gray, March 19, 2000, DFCI
c Copyright (C) 2000 Robert Gray
c Distributed under the GNU public license
c
c t-statistics
c first ng1 columns of d assumed to be group 1, other ng-ng1 assumed to be 
c group2.  Note: single precision stats
c mg=max allowed in group 1
      subroutine ptn(d,n,ng,ng1,z,dm,pzi,pzm,pzml,pzmr,nperm,ix,inml,
     $     inmr)
      real d(n,ng),z(n),zt,ztm,ztml,ztmr,dm(n),tsum
      integer pzi(n),pzm(n),pzml(n),pzmr(n),inml(n),inmr(n),inmml,inmmr
      integer mg
      parameter(mg=200)
      integer n,ng,ng1,ig1(mg),ig2(2*mg),i,j,k,ng2,nperm,isamp,ix(3)
      double precision dmm,dss
c initialize
      ng2=ng-ng1
      if (ng1.gt.mg) return      
      do 5 j=1,ng1
         ig1(j)=j
 5    continue 
c t statistics
      ztml=0
      ztmr=0
      do 10 i=1,n
c standardize rows of d
         dmm=0
         do 12 j=1,ng
            dmm=dmm+d(i,j)
 12      continue 
         dmm=dmm/ng
         dss=0
         do 13 j=1,ng
            dss=dss+(d(i,j)-dmm)**2
 13      continue 
         dss=sqrt(dss/(ng-1))
         do 14 j=1,ng
            d(i,j)=(d(i,j)-dmm)/dss
 14      continue 
c since the overall pooled variance orders the statistics in the same way
c as the within group SS pooled variance, don't need to compute variances,
c and since everything else is the same for each row, just need to compute
c sum in each group
         z(i)=tsum(ng1,d(i,1),n,ig1)
         if (z(i).gt.ztmr) then
            ztmr=z(i)
            inmmr=i
         else if (z(i).lt.ztml) then
            ztml=z(i)
            inmml=i
         endif
         inmr(i)=0
         inml(i)=0
 10   continue 
      if (nperm.le.0) then
         isamp=0
         nperm=1
         inmr(inmmr)=1
         inml(inmml)=1
         do 16 i=1,n
            pzi(i)=1
            pzml(i)=1
            pzmr(i)=1
            pzm(i)=1
 16      continue 
      else
         isamp=1
         do 17 i=1,n
            pzi(i)=0
            pzml(i)=0
            pzmr(i)=0
            pzm(i)=0
 17      continue 
      endif
 20   if (isamp.gt.0) then
         if (isamp.gt.nperm) go to 80
         isamp=isamp+1
         call rnsub(ng1,ng,ig1,ig2,ix)
      else  
         if (ig1(ng1).lt.ng) then
            ig1(ng1)=ig1(ng1)+1
         else
            do 25 j=ng1-1,1,-1
               if (ig1(j).lt.ng-ng1+j) then
                  ig1(j)=ig1(j)+1
                  do 26 k=j+1,ng1
                     ig1(k)=ig1(k-1)+1
 26               continue 
                  go to 27
               endif
 25         continue 
c no more permutations
            go to 80
         endif
 27      nperm=nperm+1
      endif
      ztm=0
      ztml=0
      ztmr=0
      do 40 i=1,n
         zt=tsum(ng1,d(i,1),n,ig1)
         if (abs(zt).ge.abs(z(i))) pzi(i)=pzi(i)+1
         if (abs(zt).ge.ztm) ztm=abs(zt)
         if (zt.lt.ztml) then
            ztml=zt
            inmml=i
         else if (zt.gt.ztmr) then
            ztmr=zt
            inmmr=i
         endif
 40   continue 
      do 45 i=1,n
         if (ztm.ge.abs(z(i))) pzm(i)=pzm(i)+1
         if (ztml.le.z(i)) pzml(i)=pzml(i)+1
         if (ztmr.ge.z(i)) pzmr(i)=pzmr(i)+1
 45   continue 
      inml(inmml)=inml(inmml)+1
      inmr(inmmr)=inmr(inmmr)+1
      go to 20
 80   continue
c calculate real stats so users don't get confused
      ng2=ng-ng1
      do 61 i=1,n
         call tst2(d(i,1),ng1,ng2,n,z(i),dm(i))
 61   continue 
      return
      end

      function tsum(ng,d,n,ig1)
      real d(n,*),tsum
      integer n,ng,ig1(ng),i
      double precision dd
      dd=0
      do 10 i=1,ng
         dd=dd+d(1,ig1(i))
 10   continue 
      tsum=dd
      return
      end

      subroutine tst2(d,ng1,ng2,n,tst,dm)
c columns 1 to ng1 in group 1, ng1+1 to ng1+ng2 in group 2
      real d(n,ng1+ng2),tst,dm
      double precision dm1,dm2,dss1,dss2
      integer ng1,ng2,n,i
      dm1=0
      dm2=0
      dss1=0
      dss2=0
      do 10 i=1,ng1
         dm1=dm1+d(1,i)
 10   continue 
      dm1=dm1/ng1
      do 11 i=1,ng1
         dss1=dss1+(d(1,i)-dm1)**2
 11   continue 
      do 12 i=1,ng2
         dm2=dm2+d(1,ng1+i)
 12   continue 
      dm2=dm2/ng2
      do 13 i=1,ng2
         dss2=dss2+(d(1,ng1+i)-dm2)**2
 13   continue 
      dm=dm1-dm2
      if (dss1.eq.0.and.dss2.eq.0) then
         tst=0
         return
      endif
c intermediate calculations in dp, so stats with many ties give same sp result
c regardless of order of calculations
      tst=(dm1-dm2)/sqrt((1.d0/ng1+1.d0/ng2)*(dss1+dss2)/(ng1+ng2-2))
      return
      end
	
c generate a random subset of ngp items from the population 1,...,ntot
c ix=seeds for Wichman Hill generator
      subroutine rnsub(ngp,ntot,index,itmp,ix)
      integer ngp,ntot,index(ngp),ix(3),itmp(ntot),i,ii
      real rand
      do 10 i=1,ntot
         itmp(i)=i
 10   continue 
      do 20 i=1,ngp
         ix(1) = mod(171 * ix(1), 30269)
         ix(2) = mod(172 * ix(2), 30307)
         ix(3) = mod(170 * ix(3), 30323)
         rand = mod(float(ix(1)) / 30269. + float(ix(2)) / 30307. +
     +        float(ix(3)) / 30323., 1.0)
         ii=rand*(ntot-i+1)+1
         if (ii.lt.1.or.ii.gt.ntot-i+1) write (6,*) 'error in rnsub'
         index(i)=itmp(ii)
         itmp(ii)=itmp(ntot-i+1)
 20   continue 
      return
      end

c counts complete separation
c first ng1 columns of d assumed to be group 1, other ng-ng1 assumed to be 
c group2.  Note: single precision and integers only
c mg=max allowed per group
      subroutine ptc(d,n,ng,ng1,ics,nperm,dtcs,ix)
      real d(n,ng)
      integer dtcs(4)
      integer mg,ix(3),nperm,it,ics(n)
      parameter(mg=200)
      integer n,ng,ng1,ig1(mg),ig2(2*mg),i,j,k,ng2,isamp,icmps
c initialize
      ng2=ng-ng1
      if (ng1.gt.mg.or.ng.gt.2*mg) return
      do 5 j=1,ng1
         ig1(j)=j
 5    continue 
      do 6 j=1,ng2
         ig2(j)=ng1+j
 6    continue 
      dtcs(1)=0
      do 10 i=1,n
         call tstatc(d(i,1),ig1,ng1,ig2,ng2,n,icmps)
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
         call rnsub(ng1,ng,ig1,ig2,ix)
      else  
         if (ig1(ng1).lt.ng) then
            ig1(ng1)=ig1(ng1)+1
         else
            do 25 j=ng1-1,1,-1
               if (ig1(j).lt.ng-ng1+j) then
                  ig1(j)=ig1(j)+1
                  do 26 k=j+1,ng1
                     ig1(k)=ig1(k-1)+1
 26               continue 
                  go to 27
               endif
 25         continue 
c no more permutations
            go to 80
         endif
 27      k=1
         do 30 j=1,ng
            if (k.le.ng1) then
               if (j.eq.ig1(k)) then
                  k=k+1
                  go to 30
               endif
            endif
            ig2(j-k+1)=j
 30      continue 
         nperm=nperm+1
      endif
      it=0
      do 40 i=1,n
         call tstatc(d(i,1),ig1,ng1,ig2,ng2,n,icmps)
         it=it+icmps
 40   continue 
c      write (6,*) it
      if (it.ge.dtcs(1)) dtcs(2)=dtcs(2)+1
      dtcs(3)=dtcs(3)+it
      if (it.gt.0) dtcs(4)=dtcs(4)+1
      go to 20
 80   return
      end

c determine whether this sample has complete separation
      subroutine tstatc(d,ig1,ng1,ig2,ng2,n,icmps)
      real d(n,ng1+ng2),di1,da1,di2,da2
      integer ig1(ng1),ng1,ig2(ng2),ng2,n,i,icmps
      di1=d(1,ig1(1))
      da1=di1
      di2=d(1,ig2(1))
      da2=di2
      do 10 i=1,ng1
         di1=min(di1,d(1,ig1(i)))
         da1=max(da1,d(1,ig1(i)))
 10   continue 
      do 12 i=1,ng2
         di2=min(di2,d(1,ig2(i)))
         da2=max(da2,d(1,ig2(i)))
 12   continue 
      icmps=0
      if (da1.lt.di2.or.di1.gt.da2) icmps=1
      return
      end

c routine for applying separate permutations to rows of a
C on input u is a long vector of uniform random numbers
      subroutine pa(a,n,nc,u)
      double precision a(n,nc),u(n*(nc-1)),tmp
      integer i,j,k,l,n,nc
      l=0
      do 10 i=1,n
         do 11 j=1,nc-1
            l=l+1
            k=int(u(l)*(nc-j+1))+j
            tmp=a(i,j)
            a(i,j)=a(i,k)
            a(i,k)=tmp
 11      continue 
 10   continue 
      return
      end
      
c correlation-statistics
c phen is a continuous phenotype
c statistics are the correlations between the rows of d and phen
c Note: single precision stats
c mg=max ng allowed
      subroutine ptcor(d,n,ng,phen,z,pzi,pzm,pzml,pzmr,nperm,ix,inml,
     $     inmr)
      real dip,d(n,ng),phen(ng),z(n),zt,ztm,ztml,ztmr
      integer pzi(n),pzm(n),pzml(n),pzmr(n),inml(n),inmr(n),inmml,inmmr
      integer mg
      parameter(mg=400)
      integer n,ng,ig(mg),i,j,k,nperm,isamp,ix(3)
      double precision dm,ds,tmp
c standardize phen 
      dm=0
      do 316 j=1,ng
         dm=dm+phen(j)
 316  continue 
      dm=dm/ng
      ds=0
      do 317 j=1,ng
         ds=ds+(phen(j)-dm)**2
 317  continue 
      ds=sqrt(ds)
      do 318 j=1,ng
         phen(j)=(phen(j)-dm)/ds
 318  continue 
c standardize rows of d and compute correlations, also initialize
      ztml=0
      ztmr=0
      do 311 i=1,n
         dm=0
         do 312 j=1,ng
            dm=dm+d(i,j)
 312     continue 
         dm=dm/ng
         ds=0
         do 313 j=1,ng
            ds=ds+(d(i,j)-dm)**2
 313     continue 
         ds=sqrt(ds)
         do 314 j=1,ng
            d(i,j)=(d(i,j)-dm)/ds
 314     continue 
         z(i)=dip(ng,d(i,1),n,phen,1)
         if (z(i).gt.ztmr) then
            ztmr=z(i)
            inmmr=i
         else if (z(i).lt.ztml) then
            ztml=z(i)
            inmml=i
         endif
         inmr(i)=0
         inml(i)=0
 311  continue
      if (nperm.le.0) then
         isamp=0
         nperm=1
         inmr(inmmr)=1
         inml(inmml)=1
         do 11 i=1,n
            pzi(i)=1
            pzml(i)=1
            pzmr(i)=1
            pzm(i)=1
 11      continue 
         do 5 j=1,ng
            ig(j)=1
 5       continue 
      else
         isamp=1
         do 12 i=1,n
            pzi(i)=0
            pzml(i)=0
            pzmr(i)=0
            pzm(i)=0
 12      continue 
      endif
 20   if (isamp.gt.0) then
         if (isamp.gt.nperm) go to 80
         isamp=isamp+1
         call rnperm(ng,phen,ix)
      else  
         do 397 j=2,ng-1
            tmp=phen(ng)
            do 396 k=ng,ng-j+2,-1
               phen(k)=phen(k-1)
 396        continue 
            phen(ng-j+1)=tmp
            if (ig(j).lt.j) then
               ig(j)=ig(j)+1
               nperm=nperm+1
               go to 395
            endif
            ig(j)=1
 397     continue 
         if (ig(ng).ge.ng) go to 80
         tmp=phen(ng)
         do 394 k=ng,2,-1
            phen(k)=phen(k-1)
 394     continue 
         phen(1)=tmp
         ig(ng)=ig(ng)+1
         nperm=nperm+1
      endif
 395  ztm=0
      ztml=0
      ztmr=0
c      write (6,*) nperm,(phen(k),k=1,ng),(ig(k),k=1,ng)
      do 40 i=1,n
         zt=dip(ng,d(i,1),n,phen,1)
         if (abs(zt).ge.abs(z(i))) pzi(i)=pzi(i)+1
         if (abs(zt).ge.ztm) ztm=abs(zt)
         if (zt.lt.ztml) then
            ztml=zt
            inmml=i
         else if (zt.gt.ztmr) then
            ztmr=zt
            inmmr=i
         endif
 40   continue 
      do 45 i=1,n
         if (ztm.ge.abs(z(i))) pzm(i)=pzm(i)+1
         if (ztml.le.z(i)) pzml(i)=pzml(i)+1
         if (ztmr.ge.z(i)) pzmr(i)=pzmr(i)+1
 45   continue 
      inml(inmml)=inml(inmml)+1
      inmr(inmmr)=inmr(inmmr)+1
      go to 20
 80   continue
      return
      end

c inner product. increment x by icx, y by icy
      function dip(ng,x,icx,y,icy)
      double precision u
      real dip,x(ng*icx),y(ng*icy)
      integer ng,icx,icy,ix,iy,i
      u=0
      ix=1
      iy=1
      do 10 i=1,ng
         u=u+dble(x(ix))*dble(y(iy))
         ix=ix+icx
         iy=iy+icy
 10   continue
      dip=u
      return
      end
	
c randomly permute u
c ix=seeds for Wichman Hill generator
      subroutine rnperm(ng,u,ix)
      integer ng,ix(3),i,ii
      real rand,u(ng),tmp
      do 20 i=1,ng-1
         ix(1) = mod(171 * ix(1), 30269)
         ix(2) = mod(172 * ix(2), 30307)
         ix(3) = mod(170 * ix(3), 30323)
         rand = mod(float(ix(1)) / 30269. + float(ix(2)) / 30307. +
     +        float(ix(3)) / 30323., 1.0)
         ii=rand*(ng-i+1)+1
         if (ii.lt.1.or.ii.gt.ng-i+1) write (6,*) 'error in rnperm'
         tmp=u(ng-i+1)
         u(ng-i+1)=u(ii)
         u(ii)=tmp
 20   continue 
      return
      end
