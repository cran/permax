c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

c ig1=vector of column indicators for group 1
c ng1=# elements in group 1
c column numbers in clurrent cluster run from i1min to ng
c iflg=1 means all possibilities in current cluster exhausted and ig1 reset
c to initial configuration
      subroutine upc(ig1,ng1,ng,i1min,iflg)
      integer ig1(ng1),ng1,ng,iflg,i1min,j,k
      iflg=0
      if (ig1(ng1).lt.ng) then
         ig1(ng1)=ig1(ng1)+1
      else
         do 25 j=ng1-1,1,-1
            if (ig1(j).lt.ng-ng1+j) then
               ig1(j)=ig1(j)+1
               do 26 k=j+1,ng1
                  ig1(k)=ig1(k-1)+1
 26            continue 
               return
            endif
 25      continue 
c no more permutations
         iflg=1
         do 10 j=1,ng1
            ig1(j)=i1min+j-1
 10      continue 
      endif
      return
      end
