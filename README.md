[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/sdm)](https://cran.r-project.org/package=sdm) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/babaknaimi/sdm?branch=master&svg=true)](https://ci.appveyor.com/project/babaknaimi/sdm) [![Travis-CI Build Status](https://travis-ci.org/babaknaimi/sdm.svg?branch=master)](https://travis-ci.org/babaknaimi/sdm)

[![](https://cranlogs.r-pkg.org/badges/sdm)](https://cran.r-project.org/package=sdm)

[![](http://cranlogs.r-pkg.org/badges/grand-total/sdm?color=orange)](https://cran.r-project.org/package=sdm)


# sdm
sdm is an object-oriented, reproducible and extensible R platform for species distribution modelling. The sdm package is designed to create a comprehensive modelling and simulation framework that: 1) provides a standardised and unified structure for handling species distributions data and modelling techniques (e.g. a unified interface is used to fit different models offered by different packages); 2) is able to support markedly different modelling approaches; 3) enables scientists to modify the existing methods, extend the framework by developing new methods or procedures, and share them to be reproduced by the other scientists; 4) handles spatial as well as temporal data for single or multiple species; 5) employs high performance computing solutions to speed up modelling and simulations, and finally; 6) uses flexible and easy-to-use GUI interface. For more information, check the published paper by Naimi and Araujo (2016) in the journal of Ecography.

## Installing sdm and all the required packages

sdm is also available on CRAN, therefore, it can simply be installed using the standard install.packages function as:

install.packages('sdm') 

The version of the package available on GitHub may be newer than the one available on CRAN, that can be installed using the install_github function from the devtools package:

devtools::install_github("babaknaimi/sdm")


Depending on the methods selected for modelling of species distribution, several packages may be required, that are needed to be installed on your machine before running the modelling procedure. A quick way to install all the required packages (as a guarantee of having access to the full functionaliy of sdm), is to simply use the function **installAll** offered by the sdm package. You can simply call it without any argument:

installAll()

## sdm website:

The official website of the package is currently available at http://biogeoinformatics.org

There is a Google group/forum for the users of the package where the questions can be posted and discussed:

https://groups.google.com/d/forum/rsdm


The authors of the package (Miguel Araujo & Babak Naimi) organise a summer school (usually for PhDs and Postdocs) every year, where both the conceptual theories and practices of species distribution modelling are well discussed. If you are interested, you can find more information on the next course at https://www.maraujolab.eu; see for example: 

https://www.maraujolab.eu/2025/02/13/2025-species-distributions-modelling-course/


### More information about the authors of the package:

Babak Naimi: https://www.biogeoinformatics.org/about-me/

Miguel Araujo: https://www.maraujolab.eu/people/miguel-araujo/


There is another website of the package developer (Babak Naimi) where you can find some articles about his research: https://r-gis.net
