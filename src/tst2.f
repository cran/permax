c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

      subroutine tst2(d,ng1,ng2,n,tst)
c columns 1 to ng1 in group 1, ng1+1 to ng1+ng2 in group 2
      real d(n,ng1+ng2),tst(n),dm
      double precision dm1,dm2,dss1,dss2
      integer ng1,ng2,n,i,j
      do 8 j=1,n
         dm1=0
         dm2=0
         dss1=0
         dss2=0
         do 10 i=1,ng1
            dm1=dm1+d(j,i)
 10      continue 
         dm1=dm1/ng1
         do 11 i=1,ng1
            dss1=dss1+(d(j,i)-dm1)**2
 11      continue 
         do 12 i=1,ng2
            dm2=dm2+d(j,ng1+i)
 12      continue 
         dm2=dm2/ng2
         do 13 i=1,ng2
            dss2=dss2+(d(j,ng1+i)-dm2)**2
 13      continue 
         dm=dm1-dm2
         if (dss1.eq.0.and.dss2.eq.0) then
            tst(j)=0
         else
            tst(j)=(dm1-dm2)/sqrt((1.d0/ng1+1.d0/ng2)*(dss1+dss2)/
     $           (ng1+ng2-2))
         endif
 8    continue 
c intermediate calculations in dp, so stats with many ties give same sp result
c regardless of order of calculations
      return
      end
