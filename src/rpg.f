c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license

      subroutine rpg(ng1,ng,nclust,mclust,mct1,ig1,ig2,ix)
      integer ng1,ng,ix(3),nclust,mclust(nclust), mct1(nclust),ig1(ng1)
      integer ig2(ng),i,j,k,l,m
      if (nclust.le.1) then
         call rnsub(ng1,ng,ig1,ig2,ix)
      else
         k=0
         l=0
         m=0
         do 320 j=1,nclust
            call rnsub(mct1(j),mclust(j),ig1(l+1),ig2(m+1),ix)
            do 321 i=1,mct1(j)
               ig1(l+i)=ig1(l+i)+k
 321        continue 
            do 322 i=1,mclust(j)-mct1(j)
               ig2(m+i)=ig2(m+i)+k
 322        continue 
            l=l+mct1(j)
            k=k+mclust(j)
            m=m+mclust(j)-mct1(j)
 320     continue 
      endif
      return
      end
