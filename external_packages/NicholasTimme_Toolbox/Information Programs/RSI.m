function [RSIOut] = RSI(CountsMat)
%RSI calculates the redundancy-synergy index.
%   [RSIOut] = RSI(CountsMat) is the redundancy-synergy index for the
%   variables in CountsMat. It is based the 2001 paper by Chechik:
%
%   G. Chechik, A. Globerson, N. Tishby, M. J. Anderson, E. D. Young, and 
%   I. Nelken, in Neaural Information Processing Systems 14, Vol. 1, 
%   edited by T. G. Dietterich, S. Becker, and Z. Ghahramani (MIT Press, 
%   2001) p. 173.
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
%   RSIOut: The redundancy-synergy index.
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

% Find the number of X variables
N = ndims(CountsMat) - 1;

% Convert the CountsMat to a joint probability distribution. (Note, this
% will have no effect if the CountsMat is already the joint probability
% distribution.)
Pxy = CountsMat/sum(CountsMat(:));

% Calculate the sum of mutual informations between the Y variable and each
% X variable

RSIOut = 0;

for iX = 1:N
    
    % Find the X variables we will sum over
    ToElim = setdiff(1:N,iX);
    
    % Create a temporary counts matrix and permute it to isolate Y and iX.
    PxyTemp = permute(Pxy,[1,iX + 1,ToElim + 1]);
    
    % Sum over the other X variables
    PxyTemp = squeeze(sum(PxyTemp(:,:,:),3));
    
    % Calculate the mutual information and add it to the sum
    RSIOut = RSIOut + MutualInfo(PxyTemp);
    
end

% Calculate the final RSI by subtracting it from the overall mutual
% information
RSIOut = MutualInfo(CountsMat) - RSIOut;


end




