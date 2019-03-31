function [VS] = VSyn(CountsMat)
%VSyn calculates Varadan's Synergy measure for the variables in CountsMat.
%   [VS] = VSyn(CountsMat) is Varadan's Synergy measure for the variables
%   in CountsMat. It utilized the MutualInfo program and it only allows up
%   to 4 X variables. It is based on the 2006 paper by Varadan:
%
%   V. Varadan, D. M. Miller III, and D. Anastassiou, Bioinformatics 22,
%   e497 (2006).
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
%   VS: Varadan's Synergy measure.
%
%
%       Version 2.0

% Version Information
%
%   1.0: 10/6/11 - The original version of the program was created before
%   and modified up to this data. (Nick Timme)
%
%   2.0: 3/21/13 - The formatting of the program was modified for inclusion
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


% Obtain the number of X variables
N = ndims(CountsMat)-1;


if N==1
    
    VS=0;
    
elseif N==2
    
    %Varadan's Synergy measure is equal to the interaction info when there
    %are 2 R variables
    VS=IntInfoMcGill(CountsMat);
    
elseif N==3
    
    % Calculate the necessary information terms
    I1=MutualInfo(squeeze(sum(sum(CountsMat,3),4)));
    I2=MutualInfo(squeeze(sum(sum(CountsMat,2),4)));
    I3=MutualInfo(squeeze(sum(sum(CountsMat,2),3)));
    I12=MutualInfo(squeeze(sum(CountsMat,4)));
    I23=MutualInfo(squeeze(sum(CountsMat,2)));
    I13=MutualInfo(squeeze(sum(CountsMat,3)));
    
    % Compute Varadan's Synergy
    VS=MutualInfo(CountsMat)-max([I1+I2+I3,I1+I23,I2+I13,I3+I12]);
    
elseif N==4
    
    % Compute the necessary information terms
    I1=MutualInfo(squeeze(sum(sum(sum(CountsMat,3),4),5)));
    I2=MutualInfo(squeeze(sum(sum(sum(CountsMat,2),4),5)));
    I3=MutualInfo(squeeze(sum(sum(sum(CountsMat,2),3),5)));
    I4=MutualInfo(squeeze(sum(sum(sum(CountsMat,2),3),4)));
    
    I12=MutualInfo(squeeze(sum(sum(CountsMat,4),5)));
    I13=MutualInfo(squeeze(sum(sum(CountsMat,3),5)));
    I14=MutualInfo(squeeze(sum(sum(CountsMat,3),4)));
    I23=MutualInfo(squeeze(sum(sum(CountsMat,2),5)));
    I24=MutualInfo(squeeze(sum(sum(CountsMat,2),4)));
    I34=MutualInfo(squeeze(sum(sum(CountsMat,2),3)));
    
    I123=MutualInfo(squeeze(sum(CountsMat,5)));
    I124=MutualInfo(squeeze(sum(CountsMat,4)));
    I134=MutualInfo(squeeze(sum(CountsMat,3)));
    I234=MutualInfo(squeeze(sum(CountsMat,2)));
    
    % Calculate Varadan's Synergy
    VS=MutualInfo(CountsMat)-max([I1+I2+I3+I4,I1+I234,I2+I134,...
        I3+I124,I4+I123,I1+I2+I34,I1+I3+I24,I1+I4+I23,I2+I3+I14,...
        I2+I4+I13,I3+I4+I12,I12+I34,I13+I24,I14+I23]);
    
elseif N>4
    
    % The script is unable to perform the calculation for any N.
    error('N is too large for VSyn');

end


end



