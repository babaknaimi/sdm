# Author: Babak Naimi, naimi.b@gmail.com
# Date (last update):  Jan 2024
# Version 1.8
# Licence GPL v3


#------
.Jackard <- function(x,y) {
  sum(apply(data.frame(x,y),1,min)) / sum(apply(data.frame(x,y),1,max))
}
#----------

