################################################################################
###                 Multivariate Information Theory Programs                 ###
################################################################################



### Introduction: ###

This toolbox contains several MATLAB programs that can be used to calculate 
various multivariate information values. 



### Authors: ###

Nicholas Timme (NMT): nmtimme@umail.iu.edu



### Variable Notation: ###

All of the programs in this toolbox are designed to calculate information 
values between different variables. There are primarily two groups of 
variables. First, there is the Y variable (often referred to as S in the 
literature). The other variables are collectively labelled as X variables (X1, 
X2, and so forth). (The X variables are often referred to as R variables in 
the literature.)



### Units: ###

Unless noted, all information values are in units of bits.



### Discrete Probability Distributions: ###

At this time, the information programs are only capable of calculating 
information values for discrete probability distributions. 



### Standard Function Inputs: ###

Inputs to the information programs are rank N+1 tensors (where N is the number 
of X variables). The elements of the tensors are either joint probability values
or counts for a joint state (though the tensor must be exclusively composed of
probability values or counts values). The first index of the tensor corresponds
to the state of the Y variable. The second through N+1 indexes of the tensor
correspond to the states of the X1 to XN variables. Consider the following 
example:

Suppoes we have a Y variable and two X variables. Let Y have 3 states 
(Y=1,2,3), X1 have 2 states (X1=1,2), and X2 have 2 states (X2=1,2).

Suppose the state Y=3, X1=1, X2=2 occurs 17 times, then element (3,1,2) of the 
tensor will be 17. 

For instance, in MATLAB, we might produce the following tensor (We have 
labelled the input tensor as 'Counts' as it contains the number of counts for 
specific states. We will follow this convention throughout.):

Counts(:,:,1) =

    18     2
    10     8
    16    19


Counts(:,:,2) =

    16     0
    20    12
    17    19

So, from this counts tensor, we can infer that the state Y=3, X1=2, X2=1 
occured 19 times. Alternatively, instead of using an input that contains 
counts for the various states, we could input the joint probability distrbution
created using the following MATLAB command:

PJoint = Counts/sum(Counts(:));

Either inputs, state counts or joint probability distributions, are acceptable
to the information programs. 



### Test Data: ###

This toolbox also contains a small portion of test data. The file 
2InputTestData.mat contains test data and the necessary matrices (see below) 
for the partial information lattice program. Specifically, it contains counts 
arrays for several simple 2 input logic gates and a random counts array for 
2 X variables and 1 Y variable.



### File List: ###

#Partial Information Decomposition#

(These programs calculate the partial information decomposition terms between 
one Y variable and a set of N X variables. N must be greater than or equal to 
2. Theoretically, N could range to infinity, but practically, the programs are 
limited to N less than 5.)

P. L. Williams and R. D. Beer, arXiv:1004.2515v1 (2010)

PIDLattice: Creates matrices that contain information about the lattice 
	used in the partial information decomposition. It needs to be run only 
	once for a given number of X variables. It takes as its input the 
	number of X variables (N).

PID: Calculates the partial information terms. It takes as its input the 
	counts array and the matrices created by the PIDLattice program.

PIDTermLabels: Shows the label for each term in the partial information 
	decomposition. It allows the user to identify the desired term/s 
	from the output of PartialInfo. It takes as its input the number 
	of X variables and it uses the Lattice program.


#Mutual Information#

(These programs calculate the mutual information between various groups of 
variables.)

MutualInfo: Calculates the mutual information between the Y variable and 
	all of the X variables considered as one vector values variable 
	{X1, ... XN}. It takes as its input the counts array. 

MutualInfoPairs: Calculates the mutual informations between the Y variable 
	and each of the X variables individually. It produces N mutual 
	information values. It takes as its input the counts array.


#Interaction Information#

(This program calculates the interaction information as described by McGill 
in his 1954 Psychometrika article.)

W. J. McGill, Psychometrika 19, 97 (1954)

IntInfoMcGill: Calcuates the interaction information between all variables. 
	It does not differentiate between Y and the other variables. It 
	utilizes the expansion for the interaction information described by 
	Jakulin 2003. It takes as its input the counts array.


#Total Correlation#

(This program calculates the total correlation as described by Watanabe in 
1960.)

S. Watanabe, IBM JJ. Res. Dev. 4, 66, (1960)

TotalCor: Calculates the total correlation between all variables. It does 
	not differentiate between Y and the other variables. It takes as 
	its input the counts array.


#Dual Total Correlation#

(This program calculates the dual total correlation as described by Watanabe
in 1978.)

T. Han, Information and Control 36, 133 (1978).

DualTotalCor: Calculates the dual total correlation between all
	variables. It does not differentiate between Y and the other variables.
	It takes as its input the counts array.


#Varadan's Synergy#

(This program calculates the synergy measure proposed by Varadan in 2006.)

V. Varadan, D. M. Miller III, and D. Anastassiou, Bioinformatics 22, e497 
(2006)

VSyn: Calculates the value of Varadan's synergy measure for all variables. 
	It does not differentiate between Y and the X variables. It can accept 
	up to four X variables. It utilizes the MutualInfo program. It takes as
	its input the counts matrix. 


#Redundancy-synergy Index#

(This program calculates the redundancy measure introduced by Chechik in 2001.)

G. Chechik, A. Globerson, N. Tishby, M. J. Anderson, E. D. Young, and I. Nelken,
in Neural Information Processing Systems 14, Vol. 1, edited by T. G. Dietterich,
S. Becker, and Z. Ghahramani (MIT Press, 2001) p. 173.

	RSI: Calculates the value of the redundancy-synergy index between 
	the Y and all of the X variables. It uses the MutualInfo program. It 
	takes as its input the counts matrix.


#Delta I#

(This program calculates the Delta I measure proposed by Latham and Nirenberg 
in 2001.)

S. Nirenberg, S. M. Carcieri, A. L. Jacobs, and P. E. Latham, Nature 411, 698 
(2001)

P. E. Latham and S. Nirenberg, Journal of Neuroscience 25, 5195 (2005)

DeltaI: Calculates the value of Delta I between the Y variable and all of the 
	X variables. It takes as its input the counts matrix.


#Entropy#

(These programs calculate various entropy quantities.)

Entropy: Calculates the entropy for the probability distributions 
	associated with the input counts matrix. It does not differentiate 
	between the X and Y variables. It takes as its input the counts matrix.

EntropyY: Calculates the entropy of the Y variable. It takes as its 
	input the counts matrix.

EntropyX: Calculates the entropy of the X variables considered as a
	single vector valued variable.

EntropyXGivY: Calculates the conditional entropy of the X 
	variables considered as a single vector valued variable given the Y
	Variable.

EntropyYGivX: Calculates the conditional entropy of the Y 
	variable given the X variables considered as a single vector valued 
	Variable.
