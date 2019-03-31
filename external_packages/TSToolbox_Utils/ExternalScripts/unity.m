% function U = unity(A) 
% A is n x m  matrix, and each column is a signal.
% the function gives a unity signal for each column
% i.e. [signal-mean(signal)] / std(signal)

function U = unity(A)

meanA = mean(A);
stdA = std(A);

U = (A - repmat(meanA,size(A,1),1))./repmat(stdA,size(A,1),1);


