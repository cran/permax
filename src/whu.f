c Copyright (C) 2000, 2002 Robert Gray
c Distributed under the GNU public license
	
      function whu(ix)
      real whu
      integer ix(3)
      ix(1) = mod(171 * ix(1), 30269)
      ix(2) = mod(172 * ix(2), 30307)
      ix(3) = mod(170 * ix(3), 30323)
      whu = mod(float(ix(1)) / 30269. + float(ix(2)) / 30307. +
     +     float(ix(3)) / 30323., 1.0)
      return
      end
