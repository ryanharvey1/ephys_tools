################################################################################

###                            Complexity Package                            ###
################################################################################




### Authors ###



Najja Marshall: njm2149@columbia.edu

Nicholas Timme: nicholas.m.timme@gmail.com




### Introduction ###



This folder contains programs that can be used to estimate the neural complexity (Tononi et al., 1994) for a system using its spike raster.  Use help or see the code for each function for additional information.

The package contains the following programs (organized by functionality):



- asdf2toraster
- complexity

- demo

complexity

Also included is a sample data set (sample_data.mat) for use with the demo files. Additional information about each program is provided below.




### asdf2toraster ###



Extracts spike times from asdf2-formatted data and converts them into a full spike raster. asdf2 is a structure-array based format used in the included AvStat package.







### complexity ###



Estimates the neural complexity for a system whose activity is contained in a spike raster.  Complexity is computed as the integrated difference between an expected linear increase in the system's integrated information with subset size and the actual integrated information for j subsets of size k (see Tononi et al., 1994 for a thorough introduction).  The program can also return the integrated information for the raster.




### democomplexity ###



Demonstrates the functionality of the package by computing the neural complexity.  The program also produces a plot of the two metrics whose difference gives the complexity measure.  The spike raster is that of real experimental data, contained in the file sample_data.mat.