c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

c inner product. increment x by icx, y by icy
c in the stratified case, have normalized within strata, so the combined
c inner product is the sum over strata of the within strata correlations
c in this case the statistic is the average of the within strata corrs.
      function dip(ng,x,icx,y,icy,nclust,mclust,istrt,wght)
      double precision u,u1
      real dip,x(ng*icx),y(ng*icy),wght(nclust)
      integer ng,icx,icy,ix,iy,i,j,k,nclust,istrt,mclust(nclust)
      u=0
      ix=1
      iy=1
      if (istrt.ne.1) then
         do 10 i=1,ng
            u=u+dble(x(ix))*dble(y(iy))
            ix=ix+icx
            iy=iy+icy
 10      continue
      else
         k=0
         do 12 j=1,nclust
            u1=0
            do 13 i=k+1,k+mclust(j)
               u1=u1+dble(x(ix))*dble(y(iy))
               ix=ix+icx
               iy=iy+icy
 13         continue
            u=u+wght(j)*u1
            k=k+mclust(j)
 12      continue 
      endif
      dip=u
      return
      end
	