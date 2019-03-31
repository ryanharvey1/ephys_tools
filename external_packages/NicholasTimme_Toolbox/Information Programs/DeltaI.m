function [dI] = DeltaI(CountsMat)
%DeltaI calculates the value of Delta I.
%   [dI] = DeltaI(CountsMat) is the value of Delta I between the Y 
%   variable and all the X variables. It uses the Delta I measure 
%   introduced by Latham and Nirenberg. 
%
%   S. Nirenberg, S. M. Carcieri, A. L. Jacobs, and P. E. Latham, Nature 
%   411, 698 (2001).
%
%   P. E. Latham and S. Nirenberg, Journal of Neuroscience 25, 5195 (2005).
%
%   Inputs
%
%   CountsMat: An array that contains the counts (or joint probability 
%   values) of the various states of the variables. The first index 
%   corresponds to the state of the Y variable. The second through N+1 
%   indexes correspond to the states of the X1 to XN variables. 
%
%   Outputs
%
%   dI: Delta I.
%
%
%       Version 2.0

% Version Information
%
%   1.0: 10/6/11 - The original version of the program was created before
%   and modified up to this data. (Nick Timme)
%
%   2.0: 3/25/13 - The formatting of the program was modified for inclusion
%   in the toolbox. (Nick Timme)
%

%==============================================================================
% Copyright (c) 2013, The Trustees of Indiana University
% All rights reserved.
% 
% Authors: Nick Timme (nmtimme@umail.iu.edu)
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
%   1. Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
% 
%   2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
% 
%   3. Neither the name of Indiana University nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==========================================================================


% Obtain the number of states for each variable
nS = size(CountsMat);

% Obtain the number of X variables
N = length(nS) - 1;

% Convert the CountsMat to a joint probability distribution. (Note, this
% will have no effect if the CountsMat is already the joint probability
% distribution.)
Pxy = CountsMat/sum(CountsMat(:));


% First, calculate the actual conditional probability distribution
PyGivx = Pxy./repmat(sum(Pxy,1),[nS(1),ones([1,N])]);
PyGivx(~isfinite(PyGivx)) = 0;


% Now, calculate the marginal probability of the Y variable
Py = repmat(sum(Pxy(:,:),2),[1,nS(2:(N + 1))]);


% Now, calculate the independent conditional probabilities of the X
% variables on the Y variable

% We can assume that we have at least one X variable to start the product
PxGivyInd = repmat(sum(Pxy(:,:,:),3),[1,1,nS(3:(N + 1))]) ./ Py;
PxGivyInd(~isfinite(PxGivyInd)) = 0;

% Now include the additional X variables
for iX = 2:N
    
    % Find the X variables we will sum over
    ToElim = setdiff(1:N,iX);
    
    % Permute the X variables to isolate Y and iX
    Temp = permute(Pxy,[1,iX + 1,ToElim + 1]);
    
    % Figure out the sizes for repmat
    SizeRep = nS;
    SizeRep(1) = 1;
    SizeRep(iX + 1) = 1;
    
    % Figure out the sizes for reshape
    SizeResh = ones([1,N + 1]);
    SizeResh(SizeRep == 1) = nS(SizeRep == 1);
    
    % Combine the results for this variable with the others in the product
    Temp2 = repmat(reshape(sum(Temp(:,:,:),3),SizeResh),SizeRep) ./ Py;
    Temp2(~isfinite(Temp2)) = 0;
    PxGivyInd = PxGivyInd .* Temp2;
    
end


% Now, calculate the independent joint probability for the X variables
PxInd = repmat(sum(PxGivyInd .* Py,1),[nS(1),ones([1,N])]);


% Now, calculate the independent conditional probability of Y given the X
% variables.
PyGivxInd = (PxGivyInd .* Py) ./ PxInd;
PyGivxInd(~isfinite(PyGivxInd)) = 0;


% Now, calculate the value of delta I
temp = Pxy .* log2(PyGivx ./ PyGivxInd);
temp(~isfinite(temp)) = 0;
dI = sum(temp(:));


    




end


