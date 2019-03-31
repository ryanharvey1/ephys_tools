function [CountsMat] = BackPropNetBinner(NetAct,Rules,StopPoint,BinNum,NumIn)
%BackPropNetBinner bins the truth tables from the back-propagation network.
%These data can then be used to calculate information values using the
%multivariate information measures from the multivariate information
%toolbox. The outputs of the network will be labelled as the Y variable.
%The inputs to the network will be lablled as the X variables. The binner
%is limited to data from back-propagation networks with only one output and
%2 to 4 inputs.

% Inputs

% NetAct: the recorded outputs from the network for each truth table 
% recording throughout development. It is a rank 3 tensor. Its dimensions 
% are the number of input states (2^NumIn) by NumOut by the total number 
% of possible truth table recordings ((NumTotEvo/NumEvoSteps)+1). All
% values after the point at which the network stops developing will be
% zero.

% Rules: the training set. The first NumIn columns contain the possible
% input values. The NumIn+1 to NumIn+NumOut columns contain the desired 
% output value. Each row is one input/output mapping. So, the AND-gate
% would be expressed as follows:

% [-1, -1, -1;
%   1, -1, -1;
%  -1,  1, -1;
%   1,  1,  1]

% Note that each distinct input state should appear in the rules only once.

% StopPoint: the number of truth table recordings 

% BinNum: the number of bins for the output (Y variable) values.

% NumIn: the number of input nodes in the network.


% Output

% CountsMat is a tensor that contains the counts of the various states of the
% variables. The first index corresponds to the state of the Y variable.
% The second through N+1 indexes correspond to the states of the X1 to XN
% variables. 



% Trim off the network results after the network stopped developing
[~,~,TotalTimeSteps]=size(NetAct);
NetAct(:,:,(StopPoint+1):TotalTimeSteps)=[];

% Find the size of the bins of the Y variables. (2 is the maximum range of
% the Y variable.)
BinSize=2/BinNum;


% Pre-allocate space for the counts matrix.
if NumIn==2
    CountsMat=zeros([BinNum,2,2,StopPoint]);
elseif NumIn==3
    CountsMat=zeros([BinNum,2,2,2,StopPoint]);
elseif NumIn==4
    CountsMat=zeros([BinNum,2,2,2,2,StopPoint]);
end


% Now bin the values of the Y variables
NetAct=NetAct+1;
NetAct=ceil(NetAct/BinSize);

% Now bin the values of the X variables
Rules(:,NumIn+1)=[]; % Trim off the output values
Rules(Rules==1)=2; % The rules are already + or - 1
Rules(Rules==-1)=1;

% Now put the binned data in the counts matrix

for i=1:StopPoint
    for j=1:(2^NumIn)
        CountsMat(NetAct(j,1,i),Rules(j,1),Rules(j,2),i)=1;
    end
end
    


end

