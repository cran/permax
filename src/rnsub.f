c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license
	
c generate a random subset of ngp items from the population 1,...,ntot
c ix=seeds for Wichman Hill generator
      subroutine rnsub(ngp,ntot,index,itmp,ix)
      integer ngp,ntot,index(ngp),ix(3),itmp(ntot),i,ii
      real whu
      do 10 i=1,ntot
         itmp(i)=i
 10   continue 
      do 20 i=1,ngp
         ii=whu(ix)*(ntot-i+1)+1
         if (ii.lt.1.or.ii.gt.ntot-i+1) write (6,*) 'error in rnsub'
         index(i)=itmp(ii)
         itmp(ii)=itmp(ntot-i+1)
 20   continue 
      return
      end
