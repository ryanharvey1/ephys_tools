function [NetAct,TotAbsError,StopPoint] = BackPropNet(NumIn,NumOut,NumHidNodes,NumHidLay,Rules,InitMax,InitMin,NumEvoSteps,NumTotEvo,LearningRate,Bias,AllowedError)
%BackPropNetToPublish creates and runs a back-propagation network designed
%to reproduce the input/output relationship specificed by the Rules matrix.
%In our case, the intput/output mapping is a truth table that defines a
%logic gate. 

% Consult the following citation for a review of the back-propagation
% network:

% P. J. Werbos, Proc. IEEE 78, 1550 (1990).

% This program was written by Nicholas Timme. It was last modified on
% 8/31/11.

% In order to prevent redundant calculates, the network only develops until
% it satisfactorally reproduces the input/output mapping. Once it reaches
% the desired mapping, it stops developing. This simplification is possible
% because the desired input/output mapping is a stable attractor and the
% network will not deviate from the desired mapping relationship once it 
% reaches that mapping. If the network cannot reach its
% target mapping, then its development ceases according to the variables
% listed below.


% Inputs

% NumIn: the number of input nodes in the network.

% NumOut: the number of output nodes in the network

% NumHidNodes: the number of nodes in each hidden layer. 

% NumHidLay: the number of hidden layers in the network. This number must
% be greater than or equal to two.

% Rules: the training set. The first NumIn columns contain the possible
% input values. The NumIn+1 to NumIn+NumOut columns contain the desired 
% output value. Each row is one input/output mapping. So, the AND-gate
% would be expressed as follows:

% [-1, -1, -1;
%   1, -1, -1;
%  -1,  1, -1;
%   1,  1,  1]

% Note that each distinct input state should appear in the rules only once

% InitMax and InitMin: the limits on the initial weights. The initial 
% weights will be set randomly. Usually this range is -1 to 1.

% NumEvoSteps: the number of development steps between each pause to record
% the truth table for the network

% NumTotEvo: the maximum total number of development time steps. Each development
% time step consists of inputing a randomly chosen input into the network,
% obtaining and output from the network, and adjusting the weights in the
% network according to the rules of the back-propagation algorithm. Make
% sure NumTotEvo is divisible by NumEvoSteps.

% LearningRate: the learning rate parameter for the back-propagation
% process

% Bias: the level of the biases in the network

% AllowedError: the total allowed error between the target outputs and the
% outputs of the network. Once the network produces a truth table with an
% error less than AllowedError, the network will stop development. 



% Outputs

% NetAct: the recorded outputs from the network for each truth table 
% recording throughout development. It is a rank 3 tensor. Its dimensions 
% are the number of input states (2^NumIn) by NumOut by the total number 
% of possible truth table recordings ((NumTotEvo/NumEvoSteps)+1). All
% values after the point at which the network stops developing will be
% zero.

% TotAbsError: the total absolute error between the outputs of the network
% and the values of the outputs in the truth table for each truth table
% recorded through development. It is a vector with a length of the number
% of possible truth table recordings ((NumTotEvo/NumEvoSteps)+1). All
% values after the point at which the network stops developing will be
% zero.

% StopPoint: the number of truth table recordings 



% Activity at each node will be calculated via a linear summation followed
% by passing the sum to a sigmoid. For now, this sigmoid will not be
% alterable. The formula for the sigmoid is:

% f(x)=2./(1+exp(-x))-1

% The formula for the derivative of the sigmoid is:

% f'(x)= 2exp(-x)./(1+exp(-x)).^2



% Re-orient the rules matrix
Rules=Rules';





% Create the hidden layer weight matrix. This matrix contains the weights
% between the hidden layers. It does not contain the weights between the
% output and the last hidden layer or the inputs and the first hidden
% layer.

HiddenWeights=(InitMax-InitMin)*rand([NumHidNodes+1,NumHidNodes+1,NumHidLay-1])+InitMin;

% Create the weights matrix between the input layer and the first hidden
% layer.

InWeights=(InitMax-InitMin)*rand([NumHidNodes+1,NumIn+1])+InitMin;

% Create the weights matrix between the output nodes and the last hidden
% layer.

OutWeights=(InitMax-InitMin)*rand([NumOut,NumHidNodes+1])+InitMin;


% Pre-allocate the network activity (include the bias)

NetAct=zeros([NumOut,2^NumIn,(NumTotEvo/NumEvoSteps)+1]);


% The input will be the complete set of all possible inputs. This set is
% represented in the rules. Create a set of inputs for these states

NetActInputs=ones([NumIn+1,2^NumIn]);   % Including the biases
NetActInputs(1:NumIn,:)=Rules(1:NumIn,:);


% Pre-allocate a vector to record the absoluate value difference between
% the target result and the actual result through development.

TotAbsError=zeros([1,(NumTotEvo/NumEvoSteps)+1]);





% Run activity on the initial network with random weights

TempAct=2./(1+exp(-InWeights*NetActInputs))-1;

TempAct((NumIn+1),:)=Bias; % reset biases

for i=1:(NumHidLay-1)
    TempAct=2./(1+exp(-HiddenWeights(:,:,i)*TempAct))-1;
    TempAct((NumIn+1),:)=Bias; % reset biases
end

NetAct(1:NumOut,:,1)=2./(1+exp(-OutWeights*TempAct))-1;


% Now record the total absolute difference between target and
% actual results. Since the inputs match the rules matrix, just use that.

TotAbsError(1)=sum(abs(NetAct(1:NumOut,:,1)-Rules((NumIn+1):(NumIn+NumOut),:)));





% Now evolve the network with pauses to record the truth tables
PauseIndexBig=2;
PauseIndexSmall=1;
i=1;
while (TotAbsError(PauseIndexBig-1)>AllowedError)&&(i<=NumTotEvo)

    % First, record the truth table if we are at a pause point.
    if PauseIndexSmall==NumEvoSteps
        
        TempAct=2./(1+exp(-InWeights*NetActInputs))-1;
        TempAct((NumIn+1),:)=Bias; % reset biases
        
        for j=1:(NumHidLay-1)
            TempAct=2./(1+exp(-HiddenWeights(:,:,j)*TempAct))-1;
            TempAct((NumIn+1),:)=Bias; % reset biases
        end

        NetAct(1:NumOut,:,PauseIndexBig)=2./(1+exp(-OutWeights*TempAct))-1;
        
        % Now record the total absolute difference between target and
        % actual results

        TotAbsError(PauseIndexBig)=sum(abs(NetAct(1:NumOut,:,PauseIndexBig)-Rules((NumIn+1):(NumIn+NumOut),:)));
        
        
        PauseIndexSmall=0;
        PauseIndexBig=PauseIndexBig+1;
    end
    
    % Now evolve the network
    
    % Calculate the activity. First, select a random rule from the rules
    % list to test.
    
    HidActTemp=zeros([NumHidNodes+1,NumHidLay]);
    
    RuleNum=randi(2^NumIn);
    



    HidActTemp(:,1)=2./(1+exp(-InWeights*[Rules(1:NumIn,RuleNum);1]))-1;
    HidActTemp(NumHidNodes+1,1)=Bias; % reset bias
    
    for j=2:NumHidLay
        HidActTemp(:,j)=2./(1+exp(-squeeze(HiddenWeights(:,:,(j-1)))*HidActTemp(:,(j-1))))-1;
        HidActTemp(NumHidNodes+1,j)=Bias; % reset bias
    end
    
    Output=2./(1+exp(-OutWeights*HidActTemp(:,NumHidLay)))-1;
    
  
    % Now calculate the error vector for the output
    
    OutputError=(Rules((NumIn+1):(NumIn+NumOut),RuleNum)-Output).*(2*exp(-Output)./((1+exp(-Output)).^2));

    
    
    % Evolve the output weights (this involves an outer product)
    OutWeightsNew=OutWeights+LearningRate*OutputError*HidActTemp(:,NumHidLay)';
    

    % Calculate the error vector for the last hidden layer
    
    LastHidError=(OutWeights'*OutputError).*(2*exp(-HidActTemp(:,NumHidLay))./((1+exp(-HidActTemp(:,NumHidLay))).^2));

    
    HiddenWeightsNew=zeros(size(HiddenWeights));
    HiddenWeightsNew(:,:,(NumHidLay-1))=HiddenWeights(:,:,(NumHidLay-1))+LearningRate*LastHidError*HidActTemp(:,(NumHidLay-1))';
    
  
    % Calculate the error vector for the rest of the hidden layers
    for j=((NumHidLay-2):-1:1)
        LastHidError=(squeeze(HiddenWeights(1:NumHidNodes,:,j+1))'*LastHidError(1:NumHidNodes)).*(2*exp(-HidActTemp(:,j+1))./((1+exp(-HidActTemp(:,j+1))).^2));

        HiddenWeightsNew(:,:,j)=HiddenWeights(:,:,j)+LearningRate*LastHidError*HidActTemp(:,j)';
    end
    
    % Calculate the error vector for the first hidden layer
    FirstHidError=(squeeze(HiddenWeights(1:NumHidNodes,:,1))'*LastHidError(1:NumHidNodes)).*(2*exp(-HidActTemp(:,1))./((1+exp(-HidActTemp(:,1))).^2));


    InWeightsNew=InWeights+LearningRate*FirstHidError*[Rules(1:NumIn,RuleNum);1]';


    % Update the weights
    OutWeights=OutWeightsNew;
    HiddenWeights=HiddenWeightsNew;
    InWeights=InWeightsNew;
    
    % Increase the indexes
    PauseIndexSmall=PauseIndexSmall+1;
    i=i+1;
end

% Record the point at which the network stopped developing.
StopPoint=PauseIndexBig-1;

% Reshape the network output activity so it matches the format of the
% rules.
NetAct=reshape(NetAct,[2^NumIn,NumOut,(NumTotEvo/NumEvoSteps)+1]);


end

