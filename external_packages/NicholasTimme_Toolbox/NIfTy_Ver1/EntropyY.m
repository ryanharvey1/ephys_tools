function [ENT] = EntropyY(CountsMat)
%EntropyY calculates the entropy in bits of all the Y variable in
%CountsMat.
%   [ENT] = EntropyY(CountsMat) is the entropy of the the Y variable in
%   CountsMat.
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
%   ENT: Entropy of the Y variable.
%
%
%       Version 2.0

% Version Information
%
%   1.0: 10/6/11 - The original version of the program was created before
%   and modified up to this data. (Nick Timme)
%
%   2.0: 3/27/13 - The formatting of the program was modified for inclusion
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

% Convert the CountsMat to a joint probability distribution. (Note, this
% will have no effect if the CountsMat is already the joint probability
% distribution.)
Pxy = CountsMat/sum(CountsMat(:));

% Find the probability distribution of the Y variable
Py = sum(Pxy(:,:),2);

% Calculate the entropy
Temp = - Py .* log2(Py);
Temp(~isfinite(Temp)) = 0;
ENT = sum(Temp(:));



end

