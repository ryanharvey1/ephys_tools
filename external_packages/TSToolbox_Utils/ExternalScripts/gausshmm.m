%function [State, hmm, decode] = gausshmm(Fet, nStates, updateOM, verbose)
% train HMM software using a Gaussian observation
% model of Fet. nStates - number of mixtures/states
% updateOM - if to update the observation model
function [State, hmm, decode] = gausshmm(Fet, varargin)

[nStates, updateOM, verbose] = DefaultArgs(varargin,{2,1,1});

if isempty(strfind(path,'hmmbox'))
    addpath('/u12/antsiro/matlab/toolboxes/HMM/hmmbox');
end

netpath='/u12/antsiro/matlab/toolboxes/netlab';
if isempty(strfind(which('kmeans'),'netpath'))
	addpath(netpath);
end

[nT nDim] = size(Fet);

% Train up GMM on this data
hmm.K=nStates;

hmm=hmminit(Fet,hmm,'full');


% Train up HMM on observation sequence data using Baum-Welch
% This uses the forward-backward method as a sub-routine

hmm.train.cyc=30;
hmm.obsmodel='Gauss';
hmm.train.obsupdate=ones(1,hmm.K)*updateOM;    % update observation models ?
hmm.train.init=1;         % Yes, we've already done initialisation

hmm=hmmtrain(Fet,nT,hmm);

[decode]=hmmdecode(Fet,nT,hmm);

% Find most likely hidden state sequence using Viterbi method
State = decode(1).q_star;

