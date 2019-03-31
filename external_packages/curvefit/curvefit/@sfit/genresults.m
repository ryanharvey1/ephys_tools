function resultstrings = genresults(fitresult,goodness,output,warnstr,errstr,convmsg,clev)
%GENRESULTS   Generate cell array of strings for results window.
%
%   RESULTSSTRINGS = GENRESULTS(FITRESULT,GOODNESS,OUTPUT) takes output
%   from FIT command and pretties it up and puts it in a cell array of
%   strings.

%   Copyright 2008 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2008/12/01 07:07:12 $

if nargin<7, clev = 0.95; end

resultstrings = {};
if isempty(errstr)
    if output.exitflag < 0 
        resultstrings{end+1} = xlate('Fit computation did not converge:');
        resultstrings{end+1} = convmsg;
        resultstrings{end+1} = '';
        resultstrings{end+1} = xlate('Fit found when optimization terminated:');
        resultstrings{end+1} = '';
    elseif output.exitflag == 0
        resultstrings{end+1} = xlate('Fit computation did not converge:');
        resultstrings{end+1} = convmsg;
        resultstrings{end+1} = '';
        resultstrings{end+1} = xlate('Fit found when optimization terminated:');
        resultstrings{end+1} = '';
    end
    
    [ignore,line2,line3,line4] = makedisplay(fitresult,'f',output,clev);
    resultstrings(end+1:end+3) = {line2,line3,line4};
    
    resultstrings{end+1} = sprintf('Goodness of fit:');
    resultstrings{end+1} = sprintf('  SSE: %0.4g',goodness.sse);
    resultstrings{end+1} = sprintf('  R-square: %0.4g',goodness.rsquare);
    resultstrings{end+1} = sprintf('  Adjusted R-square: %0.4g',goodness.adjrsquare);
    resultstrings{end+1} = sprintf('  RMSE: %0.4g',goodness.rmse);

    if goodness.rsquare<0
        resultstrings{end+1} = sprintf( [
            '\n', ...
            '  Warning: A negative R-square is possible if the model\n', ...
            '           does not contain a constant term and the fit\n', ...
            '           is poor (worse than just fitting the mean).\n', ...
            '           Try changing the model or using a different StartPoint.\n'
            ] );
    end
end

