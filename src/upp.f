c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

c update permutation
c logic is convoluted, and should need less swapping, but it seems to work
c at each step permutes back to the original order before proceding to the 
c next level
      subroutine upp(phen,ng,ig,iflg)
      real phen(ng),tmp
      integer ng,ig(ng),j,k,iflg
      iflg=0
      do 397 j=2,ng-1
         tmp=phen(ng)
         do 396 k=ng,ng-j+2,-1
            phen(k)=phen(k-1)
 396     continue 
         phen(ng-j+1)=tmp
         if (ig(j).lt.j) then
            ig(j)=ig(j)+1
            return
         endif
         ig(j)=1
 397  continue 
      if (ig(ng).ge.ng) then
         iflg=1
         ig(ng)=1
      else
         ig(ng)=ig(ng)+1
      endif
      tmp=phen(ng)
      do 394 k=ng,2,-1
         phen(k)=phen(k-1)
 394  continue 
      phen(1)=tmp
      return
      end
