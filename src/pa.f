c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

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
