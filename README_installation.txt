# RChemMass is a collection of useful functions, some of which are in development
# or accessible in non-standard sources (mixed CRAN, BioConductor and github), 
# which means the dependency installation does not always work as desired.
# Herewith a short guide. 
# If you save this as .R instead of .txt, this will become a script for installing.
 
# In R
# Install Rtools:
# https://cran.r-project.org/bin/windows/Rtools/
# Be careful with the version; check the option "edit system PATH" to enable access.
# Alternative suggested by Anjana Elapavalore: 
install.packages("installr")
library("installr")
install.Rtools()

# Once Rtools is installed:
install.packages("devtools")
# this allows direct installation from e.g. github.
library(devtools)

# Ensure you have latest versions of all packages from CRAN:
# either via packages window in Rstudio or directly:
install.packages(c("curl","rsvg","enviPat","rJava", "fingerprint", "png", "rcdk", "rcdklibs"), dependencies=TRUE)

# Ensure you have the latest versions of packages available from BioConductor:
# previous Bioc install was via biocLite, now use BiocManager. 
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("MSnbase", "mzR", "RMassBank"))
# for more information see https://bioconductor.org/install/
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
# then install specific packages with
BiocManager::install(c("MSnbase", "mzR", "RMassBank"))

# if you get errors during installation about missing packages, check whether these 
# are on CRAN or BioConductor, install these seperately, test with library(package) 
# then try to install the package that failed.

# If Java issues occur (quite common with rJava and essential for rcdk functionality)
# Install the latest Java versions: https://www.java.com/en/ 
# Java 1.8 is needed for rcdk
# Ensure 32 bit and 64 bit versions are available for 64 bit systems
# Check in C:\Program Files\Java and C:\Program Files (x86)\Java 
# Test to see that library(rJava) works
library(rJava)
# If you still get error messages in R and have both Javas, try:
if (Sys.getenv("JAVA_HOME")!="")
  Sys.setenv(JAVA_HOME="")
library(rJava)
# more notes available at: https://github.com/CDK-R/cdkr 

# Test progress: 
library(RMassBank)
# if this works, the majority of the above issues should have been solved and 
# RChemMass should work (they share most of the same dependencies).

# Install the remaining dependencies from github directly:
install_github("CDK-R/rinchi")
# this is useful to create InChIs, can never be part of rcdk as it changes dirs in the background.
# to get the latest rcdk and rcdklibs 
# [WARNING - may be incompatible with some functions in current version]
library(devtools)
install_github("CDK-R/cdkr", subdir="rcdklibs")
install_github("CDK-R/cdkr", subdir="rcdk")

# External Resources:
# Some functions rely on OpenBabel, download here:
# http://openbabel.org/wiki/Main_Page (green download arrow).
# Note the path (or copy babel.exe and obabel.exe to a path of your choice).
# Functions needing this will ask for the path.

# Then finally, install RChemMass:
install_github("schymane/RChemMass")

