# KS2wrapper
This wrapper is a set of functions to collect neuroscope .xml information to use kilosort 2 in the binary files (.dat). 

If the information in the .xml file is correct, the user only needs to run the function KS2Wrapper to spike sort using kilosort 2.

If you want kilosort to not use some channels during the sorting (i.e. broken shank or dead channels), select them on neuroscope and click on skip.

This code should be self contained besides npy-matlab and Kilosort2 repositories.

## Dependencies 
- [npy-matlab](https://github.com/kwikteam/npy-matlab)
- [Kilosort2](https://github.com/MouseLand/Kilosort2)

## Installation
Download and add the folder to your MATLAB path, along with Kilosort2 and npy-matlab.

This repository is based on the [KilosortWrapper](https://github.com/brendonw1/KilosortWrapper).
