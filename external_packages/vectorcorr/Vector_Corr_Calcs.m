function RtrnCode = Vector_Corr_Calcs(AzData, DegRad, NIndepX, IndepX, ...
    IndepName, IndepAz, alfa, fidO, Plots, DataTtl, NTrials)

%Vector_Corr_Calcs.m

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Controls all the calculations and outputs
%   Set up data arrays
%   Calculate stats on individual variables
%   Compute linear-circular and circular-circular correlations
%   Generate crossplots of input data 

% Input variables:
%   AzData: input array containing variable(s)
%   DegRad: input azimuths in degrees(0) or radians(1)?
%   NIndepX: number of variables
%   IndepX: list of which columns in AzData are used as variables
%   IndepName: identifiers of variables
%   IndepAz: specify which variables are azimuths (=1)
%   alfa: significance level for tests
%   fidO: output file for summary of job (>0 implies write file)
%   Plots: specify if plots are to be made (=1)
%   DataTtl: job title
%   NTrials: number of resampling trials
% Output variable:
%   RtrnCode: error return code
%    0 = OK

% Functions and scripts called from this module:
%   Vector_Corr_SummaryStats.m
%   Vector_Corr_Analyses.m

RtrnCode = 0;
ptskip = sprintf('\n');
[NTot, nc] = size(AzData);

% Set up data files into arrays

jjj = 0;
for iii = 1:12
    if IndepX(iii) > 0
        jjj = jjj + 1;
        if IndepAz(iii) > 0 & DegRad > 0
            DataArray(:, jjj) = AzData(:, IndepX(iii))*57.3;
        else
            DataArray(:, jjj) = AzData(:, IndepX(iii));
        end
    end
end
                                              
% Calculate each variable's statistics - ThetaHat, Rsquared, KappaHat
%                                      - mean, variance

RtrnCode = Vector_Corr_SummaryStats(DataArray, fidO, IndepX, ...
                       IndepName, IndepAz);
              
% Compute correlation analyses

RtrnCode = Vector_Corr_Analyses(DataArray, IndepX, IndepName, IndepAz, ...
          alfa, fidO, DataTtl, Plots, NTrials);
if RtrnCode > 0                              
    ptstop=sprintf('JOB TERMINATING: ReturnCode = %.0f \n',RtrnCode);
    disp(ptskip); disp(ptstop); 
    if fidO > 0
        fprintf(fidO, ptskip); fprintf(fidO, ptstop);
    end
    return
end


