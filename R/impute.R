# Author: Babak Naimi, naimi.b@gmail.com
# Date of last update :  Feb 2024
# Version 1.1
# Licence GPL v3
#--------------


# Get the rarandomly selected Non-NA neighbour cell of a cell in a Raster object (p)
.getNeighCellRandom <- function(cell,p) {
  aj <- adjacent(p,cell,direction='queen',pairs=TRUE)[,2]
  x <- p[aj]
  w <- which(apply(x,1,function(x) all(!is.na(x))))
  if (length(w) > 0) {
    if (length(w) == 1) aj
    else aj[sample(w,1)]
  } else NA
}
#-------------

