c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license
    
c correlation-statistics
c phen is a continuous phenotype
c statistics are the correlations between the rows of d and phen
c Note: single precision stats
      subroutine ptcor(d,n,ng,phen,z,pzil,pzir,pzml,pzmr,nperm,ix,
     $     ig,nclust,mclust,istrt,wght,iflag,ipc,uv,zt)
      real dip,d(n,ng),phen(ng),z(n),ztml,ztmr,wght(nclust)
      real uv(nclust),zt(n)
      integer pzil(n),pzir(n),pzml(n),pzmr(n)
      integer n,ng,ig(2*ng),i,j,k,nperm,isamp,ix(3),iflag(nclust)
      integer nclust,mclust(nclust),istrt,ipc
c initialize
      if (nperm.le.0) then
         isamp=0
         nperm=1
         do 11 i=1,n
            pzil(i)=1
            pzml(i)=1
            pzmr(i)=1
            pzir(i)=1
 11      continue 
         do 5 j=1,ng
            ig(j)=1
 5       continue 
      else
         isamp=1
         do 12 i=1,n
            pzil(i)=0
            pzml(i)=0
            pzmr(i)=0
            pzir(i)=0
 12      continue 
      endif
 20   if (isamp.gt.0) then
         if (isamp.gt.nperm) go to 80
         isamp=isamp+1
         if (ipc.ne.1) then
            call rnperm(ng,phen,ix,nclust,mclust)
         else
            call rnperm(nclust,uv,ix,1,nclust)
            k=0
            do 420 j=1,nclust
               do 421 i=k+1,k+mclust(j)
                  phen(i)=uv(j)
 421           continue 
               k=k+mclust(j)
 420        continue 
            call stdms(1,phen,ng)
         endif
      else  
         if (ipc.ne.1) then
            k=1
            do 410 j=1,nclust
               call upp(phen(k),mclust(j),ig(k),iflag(j))
               if (iflag(j).lt.1) then
                  do 411 i=1,j-1
                     iflag(i)=0
 411              continue 
                  go to 27
               endif
               k=k+mclust(j)
 410        continue 
            go to 80
         else
            call upp(uv,nclust,ig,iflag)
            if (iflag(1).eq.1) go to 80
            k=0
            do 422 j=1,nclust
               do 423 i=k+1,k+mclust(j)
                  phen(i)=uv(j)
 423           continue 
               k=k+mclust(j)
 422        continue 
            call stdms(1,phen,ng)
         endif
 27      nperm=nperm+1
      endif
 395  continue
      ztmr=z(1)-1
      do 40 i=1,n
         zt(i)=dip(ng,d(i,1),n,phen,1,nclust,mclust,istrt,wght)
         if (zt(i).ge.z(i)) pzir(i)=pzir(i)+1
         if (zt(i).le.z(i)) pzil(i)=pzil(i)+1
         ztmr=max(zt(i),ztmr)
         if (ztmr.ge.z(i)) pzmr(i)=pzmr(i)+1
 40   continue 
      ztml=z(n)+1
      do 43 i=n,1,-1
         ztml=min(zt(i),ztml)
         if (ztml.le.z(i)) pzml(i)=pzml(i)+1
 43   continue 
      go to 20
 80   continue
c enforce monotonicity
      do 45 i=2,n
         pzml(i)=max(pzml(i-1),pzml(i))
         pzmr(n-i+1)=max(pzmr(n-i+1),pzmr(n-i+2))
 45   continue 
      return
      end
