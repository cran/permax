# By R Gray, June24, 2000, DFCI
# Copyright (C) 2000 Robert Gray
# Distributed under the GNU Public License (see the file COPYING)

## edit as appropriate
.First.lib <- function(a,b) {
   library.dynam('permax.so',b,a) 
}
