mle_rhythmicity
================
Welcome to mle_rhythmicity!

mle_rhythmicity is a set of MATLAB tools for analyzing the rhythmicity of 
event times. It was specifically developed for the analysis of theta 
(10 Hz) rhythmic neurons, where the most common existing methods rely on 
the binned event-time autocorrelogram. There are a number of problems with 
this technique (Climer et. al., 2015), and to overcome them we developed a 
parametric conditional-intensity function for the lags in the 
autocorrelogram window. This is considerably less biased than existing 
techniques, and allows us to do more rigorous statistics such as parameter 
estimation using the maximum-likelihood approach and to examine if 
features of rhythmicity are modulated by other covariates 
(Hinman et. al., in press).

The main two functions are mle_rhythmicity (which estimates rhythmicity 
parameters when we assume the average underlying rhythmicity is constant)
and rhythmicity_covar (which estimates rhythmicity parameters when some are
allowed to shift with a covariate). For more details, please see the 
documentation in the matlab files (doc mle_rhythmicity).

Please raise issues using the issue mechanic in the GitHub repo, or contact 
Dr. Jason Climer via email at jason.r.climer@gmail.com. 

This is a major update to the mle_rhythmicity toolset. If you'd like to get 
back to the latest version of the old code it is in the last committed 
version of the previous revision is the SHA starting with 93862ac...

Copyright (c) 2015-2016 Trustees of Boston University
All rights reserved
