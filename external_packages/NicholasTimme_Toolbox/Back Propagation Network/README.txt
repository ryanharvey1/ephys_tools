Thank you for downloading the back-propagation network toolbox. This toolbox is primarily designed to provide data to be 
analyzed by the multivariate information toolbox. Please contact me at nmtimme (at) umail.iu.edu to report bugs or offer 
suggestions on ways to improve the programs.

This toolbox contains the following programs:

BackPropNet: This program runs a back-propagation network that attempts to match the input/output mapping specified as an
input to the program. The network was designed to reproduce the input/output mappings specified by logic gates, so the 
network's development is halted periodically and the entire truth table for the network is recorded. 

The structure and developmental algorithm of the back-propagation network are described here:

P. J. Werbos, Proc. IEEE 78, 1550 (1990)

This program is able to produce mappings between any number of input and output variables. The total error between the 
output values of the network and the desired outputs is recorded also. The program is constructed to stop development
of the network once the error between the actual output values and the desired output values falls below a certain 
threshold. This decreases processing time.



BackPropNetBinner: This program takes the data from the back-propagation network and bins the values of the input and output
variables. This process is primarily intended for preparing the data for information theoretic analysis. 