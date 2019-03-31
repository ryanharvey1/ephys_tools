function [cs1,cm,cs2] = RMatrixCorrCoeffs(SS1,SB,SS2,sameTT,Stime_binsize,Btime_binsize,exclude_sameTT,rescale_timebins)
% computes the corrlation coeffients (the so-called R-matrices) of each valid cell pair of a  sleep1,behavior and sleep2 subsession triplet 
%
% part of the Explained Variance Sleep analysis script
% Version 2.2
%
% An explained variance and correlation matrix producing script 
% that works on rat data and Monkey data
% adapted by PL from Jason G

% Required:
% 1)  S-matrices for sleep 1, behavior and sleep 2: SS1, SB, SS2 
% 2)  expects a same-Tetrode flag matrix sameTT
% 3)  Be sure to set the variable time_binsize to the desired number of msec
% 4)  Set the rescale variable on or off as desired.  Rescale means that the
%     program will compute the sparsity ratio between maze and sleep matrices
%     and rescale the binsizes of sleep accordingly (idea is to account for 
%     the possibility of time compression during sleep).
% 5)  Exclude_sameTT is default to 1; this removes all within tetrode correlations
%
% Stime_binsize = 100     % binsize of sleep Q matrices in msec
% Btime_binsize = 100    % binsize of  behavior Q matrix in msec
% rescale_timebins = 0  % turn on(1)/off(0) the rescaling of timebins for S1 and S2 matrices
% exclude_sameTT = 1    % exclude cell pairs from same tetrode  



% make Q matrices (tsd) for sleep1, behaviour, sleep2

QS1 = MakeQfromS(SS1,10*Stime_binsize,'ProgressBar','none');   
QB  = MakeQfromS(SB ,10*Btime_binsize,'ProgressBar','none');
QS2 = MakeQfromS(SS2,10*Stime_binsize,'ProgressBar','none');

%  rescale time bins for S1 and S2 matrix so that their sparsities corresponds to M1
if rescale_timebins
    [AB, ABav, ABstd]  = StateVectorSparsity(QB); %calculates timescale by the sparsity in
    [AS1, AS1av, AS1std]  = StateVectorSparsity(QS1); % the matrices of sleep and maze
    [AS2, AS2av, AS2std]  = StateVectorSparsity(QS2);
   % S2_timescale_factor = 1/10   To set maze vs. sleep timescale to 1/10 as reported by Bill
   % S1_timescale_factor = 1/10
    QS1 = MakeQfromS(S_s1,10*time_binsize*S1_timescale_factor);
    QS2 = MakeQfromS(S_s2,10*time_binsize*S2_timescale_factor);
    clear AB,ABav,ABstd,AS1,AS1av,AS1std,AS2,AS2av,AS2std;
end

% make R matrices
RS1 = corrcoef(Data(QS1));
RB = corrcoef(Data(QB));
RS2 = corrcoef(Data(QS2));

%set the NaNs in the behavior data to zeros

RB(isnan(RB)) = 0;
RS1(isnan(RS1)) = 0;
RS2(isnan(RS2)) = 0;
% set all values from within tetrode comparisons to "Inf" or "-Inf"
% Note: this output will be a NaN for correlations using cells without
% spikes  
warning off;
RS1 = RS1 ./sameTT;
RB = RB ./sameTT;
RS2 = RS2 ./sameTT;
warning on;

% partial.m  these are lists of correlations% histograms of the raw cell-pair correlations

cs1  =[]; % these will contain the lower diag of the R matrices
cm  = [];
cs2 = [];

% extract the diagonals of the R matrices
N = length(SS1);   % number of spike trains in the analysis

% Creates lists of the cell-pair correlations for each epoch by
%  getting them out of the R matrix and creating a vector
for k = 1:N-1
    cs1 = full([cs1;diag(RS1,k)]);
end

for k = 1:N-1
    cm = full([cm;diag(RB,k)]);
end

for k = 1:N-1
    cs2 = full([cs2;diag(RS2,k)]);
end

% Remove any NaNs from data, this should only occur in the behavior
%  Data when the NaN from the cell not firing was replaced with a 0
%  and then the ./sameTT divides by 0, creating another NaN in the maze
%  behavior.  There really shouldn't be any NaNs in the sleep epochs

cs2(isnan(cs2)) = [];
cs1(isnan(cs1)) = [];
cm(isnan(cm)) = [];

% Remove all the "Inf" values that are "theoretically" only created by the
% "./sameTT" function done above to create an "Inf" value for within tetrode
%  comparisons

cs2(isinf(cs2)) = [];
cs1(isinf(cs1)) = [];
cm(isinf(cm)) = [];



