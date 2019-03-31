################################################################################
###                  Neuronal Avalanche Statistics Package                   ###
################################################################################




### Authors ###



Najja Marshall: njm2149@columbia.edu
Nicholas Timme: nicholas.m.timme@gmail.com

### Introduction ###



This folder contains programs that can be used to fit, analyze, and 
visualize power-law distributed data.  Additional programs in the software 
package are tailored to computing neuronal avalanche statistics.  While 
the package is largely representative of original work, the development 
of some of the all-purpose programs was inspired and guided by Aaron 
Clauset's Power Law package.  Additionally, one included program 
(rldecode) can be found freely available on the web.

 Additional information about each function can be found using help.

The package contains the following programs (organized and grouped by functionality):

- pldist

- plmle

- pvcalc

- mymnrnd

- rldecode
- gendata
- plparams

- demotruncdist


- plplottool
- demoplotting


- rastertoasdf2

- randomizeasdf2
- rebin

- avprops

- avpropvals

- sizegivdurwls

- avgshapes

- avshapecollapse

- avshapecollapsestd

- bethelattice
- brestimate
- cbmodel
- demoempdata




Also included is a sample data set (sample_data.mat) for use with one of 
the demo files. Additional information about each program is provided below.




### pldist ###



Generates perfectly power-law distributed data.  By default the program 
makes non-truncated data, but the user can include an upper and/or lower cutoff.




### plmle ###



Estimates the scaling exponent for data that is assumed to be power-law 
distributed by the method of maximum likelihood.  The program is capable 
of fitting to a truncated region given the upper and/or lower cutoff(s).




### pvcalc ###



Calculates the p-vale for a power law fit by Monte Carlo.  For 
computational efficiency, the program updates the likelihood of successful 
results using the binomial distribution and halts for statistically 
unlikely results (this feature can be turned off to enable a full 
computation).  The probability used for the binomial distribution is 
referred to as the critical p-value.

### mymnrnd ###



Generates a set of multinomially distributed random numbers.  This 
program is a custom-made version of MATLAB's mnrnd that is used for 
generating simulated data sets during the Monte Carlo process in 
pvcalc.m (above).

### rldecode ###

Run-length decoding of run-length encode data. Essentially, this program performs the inverse operation of hist. This program was written by Peter John Acklam and is freely available online.

### gendata ###

Generates integer data drawn from a collection of specified probability distributions using mymnrnd.




### plparams ###



Automatically searches for the optimal triplet of scaling exponent, 
lower and upper cutoffs for data that is assumed to be power-law 
distributed.  The program implements a refined greedy search algorithm 
to settle on a support pair.  This program is a high level macro that is 
also capable of returning the uncertainty on a fit exponent, the p-value, 
the critical p-value, and the Kolmogorov Statistic -- outputs from the 
aforementioned programs.




### demotruncdist ###



Demonstrates the capability of the software package to fit power-law 
distributed data with or without an upper and/or lower cutoff.  The 
program fits simulated data generated with pldist and also shows the 
consequences of not accounting for the cutoff values when computing the 
p-value for a distribution in order to determine whether or not it is 
power-law distributed.




### plplottool ###



Visualizes power-law distributed data.  The program is equipped to 
produce plots of data pdfs.  The program can also produce multiple, overlaid plots 
when given data in cellular arrays. The program can also overlay power-law fits.




### demoplotting ###



Demonstrates the versatility of plplot by plotting simulated power law 
data using each unique plot option.  The demo also produces two overlaid 
plots of simulated continuous data.




### rastertoasdf2 ###



Converts a spike raster into asdf2 format.  asdf2 is a structure-array 
based format that is required for avalanche statistics analysis. 

### randomizeasdf2 ###

Randomizes asdf2 format rasters using a variety of specified methods.




### rebin ###



Rebins asdf2-formatted data. 




### avprops ###



Computes the avalanche properties for a set of spike data.  Specifically, 
the program calculates the duration (number of active time bins), size 
(total number of activated neurons), shape (number of neurons activated 
at each time bin), and branching ratio (number of neurons active at time 
step t divided by the number active at time step t-1) for each avalanche.  
The program also computes the fingerprint, defined as the event times and 
frequency of occurrence of avalanches of a particular size for the data 
set.  All properties are grouped into a structure array.




### avpropvals ###



Returns the power-law values (scaling exponent, support pair, p-value) 
for the distribution of a specified avalanche property, denoted by a 
string input.  The program can analyze one of three distributions: size, 
duration, and average size given duration. Size and duration 
distributions are analyzed by means of MLE using plparams.m (above), 
while the average size given duration distribution is analyzed by means 
of weighted least squares using sizegivdurwls.m (below).




### sizegivdurwls ###



Computes the average size given duration distribution for avalanche data 
and determines its associated scaling parameter and standard deviation by 
the method of weighted least squares.




### avgshapes ###



Computes the mean temporal profiles of a given set of avalanche shapes 
and durations.  The program averages across profiles with the same 
duration.  By default the program considers all durations but an 
alternative sampling method may be indicated by the user as a variable 
argument. 




### avshapecollapse ###



Performs an automated avalanche shape collapse using average temporal 
profiles returned from avgshapes.m (above) in order to determine the 
scaling exponent 1/(sigma nu z).  The optimal exponent is determined by 
minimizing a cost function.




### avshapecollapsestd ###



Computes the uncertainty on the shape collapse estimate of 1/(sigma nu z) 
by non-parametric bootstrap.

### bethelattice ###

Generates data from a Bethe Lattice (otherwise known as a branching process). These model data can be used to test and demonstrate the analysis software, as well as in dedicated analyses.

### brestimate ###

Estimates the branching ratio for asdf2 format data under the assumption the data is sub-sampled using the method established by Wilting and Priesemann. This method allows for a better estimate of the branching ratio than simply calculating the ratio of activity at succeeding time steps (as is done in avprops.m), though it is more computationally intensive and requires additional parameters.

### cbmodel ###

Generates data from a cortical branching model. The model consists of nodes arranged on a 2-D square lattice with nearest neighbor interactions.




### demoempdata ###



Demonstrates the capability of the package to extract and analyze 
neuronal avalanches from a set of real experimental data 
(sample_data.mat).  In particular, the program produces a histogram of 
the computed avalanche branching ratios; computes the power-law 
parameters for each distribution (size, duration, and average size given 
duration); performs avalanche shape collapse and plots the results.


