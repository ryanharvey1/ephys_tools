function RtrnCode = Vector_Corr_Analyses(DataArray, IndepX, IndepName,...
           IndepAz, alfa, fidO, DataTtl, Plots, NTrials)

% Vector_Corr_Analyses.m  

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Loops over the pairs of specified variables.  For each, it 
% determines if azimuth-azimuth or linear-azimuth or linear-linear pair, 
% and calculates correlations. 
% If requested, also generates crossplots of azimuth-involved pairs.

% Input variables:
%   DataArray: input array containing azimuth (degr) and linear vars.
%   IndepX: list of which columns in AzData are used 
%   IndepName: identifiers of independent variables
%   IndepAz: indicates if variable is linear(=0) or azimuth (=1)
%   alfa: significance level for tests 
%     should be 0.10, 0.05, or 0.01 for CAssociation test
%     should be 0.10, 0.05, 0.025, 0.01, 0.005 for CircularRank test
%   fidO: output file for summary of job (>0 implies write file)   
%   DataTtl: job title
%   Plots: are plots to be generated (=1) 
%   NTrials: number of resampling trials
% Output variables:
%   RtrnCode: error return code
%    0 = OK

% Functions called from this module:
%   Vector_Corr_CircCirc1.m
%   Vector_Corr_CircTlin1.m
%   Vector_Corr_CircTmono.m
%   Vector_Corr_PMCorr.m
%   Vector_Corr_LinCirc1.m
%   Vector_Corr_CAssoc1.m
%   Vector_Corr_Plots.m

RtrnCode = 0;
                                   
ptskip = sprintf('\n');                        
pteq = sprintf('==============================================\n');
ptlin = sprintf('----------------------------------------------\n');
pt1A = sprintf('Circular-circular correlation\n');
pt1B = sprintf('Linear-circular correlation\n');
pt1C = sprintf('Linear-linear correlation\n');
disp(ptskip); disp(pteq); disp(ptskip)
if fidO > 0
    fprintf(fidO, ptskip); fprintf(fidO, pteq); fprintf(fidO, ptskip);
end

figCnt = 0;
CalcCnt = 0;
[NTot, nc] = size(DataArray);
AV1 = zeros(NTot,1); AV2 = zeros(NTot,1); LV = zeros(NTot,1);

% loop over all combinations of variables, and calculate for those
% with either or both of pair is an azimuth variable

jjj1 = 0;
for iii1 = 1:11
    if IndepX(iii1) > 0
        jjj1 = jjj1 + 1;
        
        jjj2 = jjj1;
        for iii2 = iii1+1:12
            if IndepX(iii2) > 0
                jjj2 = jjj2 + 1;
                
                % test if circular-circular correlation
                
                if (IndepAz(iii1) > 0) & (IndepAz(iii2) > 0)
                    CalcCnt = CalcCnt + 1;
                    if CalcCnt > 1
                        disp(ptlin); disp(ptskip)
                        if fidO > 0
                            fprintf(fidO, ptlin); fprintf(fidO, ptskip);
                        end
                    end
                    xlbl = IndepName(iii1,:);
                    ylbl = IndepName(iii2,:);
                    pt2 = sprintf('X%0.f  %s      X%.0f  %s \n', ...
                           iii1, xlbl, iii2, ylbl);
                    disp(pt1A); disp(pt2); disp(ptskip)
                    if fidO > 0 
                        fprintf(fidO, pt1A); fprintf(fidO, pt2);
                        fprintf(fidO, ptskip);
                    end
                    AV1 = DataArray(:,jjj1);
                    AV2 = DataArray(:,jjj2);
                    RtrnCode = Vector_Corr_CircCirc1(AV1, AV2, alfa, fidO, ...
                                                     NTrials);
                    RtrnCode = Vector_Corr_CircTmono(AV1, AV2, alfa, fidO, ...
                                                     NTrials);
                    RtrnCode = Vector_Corr_CircTlin1(AV1, AV2, alfa, ...
                                                     fidO, 0, NTrials);
                    if Plots > 0
                       xlbl = IndepName(iii1,:);
                       ylbl = IndepName(iii2,:);
                       [RtrnCode, figCntX] = Vector_Corr_Plots(1, ...
                           AV1, AV2, xlbl, ylbl, DataTtl, figCnt);
                       figCnt = figCntX;
                    end
                end
                
                % test if circular-linear / linear-circular or
                % linear-linear correlation
                
                dothis=0;
                if (IndepAz(iii1) > 0) & (IndepAz(iii2) == 0)
                    AV1 = DataArray(:,jjj1);
                    LV = DataArray(:,jjj2);
                    xlbl = IndepName(iii2,:);
                    ylbl = IndepName(iii1,:);
                    xxlbl = ylbl;  yylbl = xlbl;
                    pt3 = pt1B;
                    dothis=1;
                end
                if (IndepAz(iii1) == 0) & (IndepAz(iii2) > 0)
                    AV1 = DataArray(:,jjj2);
                    LV = DataArray(:,jjj1);  
                    xlbl = IndepName(iii1,:);  xxlbl = xlbl;
                    ylbl = IndepName(iii2,:);  yylbl = ylbl;
                    pt3 = pt1B;
                    dothis=2;
                end
                if (IndepAz(iii1) == 0) & (IndepAz(iii2) == 0)
                    LV1 = DataArray(:,jjj2);
                    LV = DataArray(:,jjj1);  
                    xxlbl = IndepName(iii1,:); 
                    yylbl = IndepName(iii2,:); 
                    pt3 = pt1C;
                    dothis=3;
                end
                
                if dothis > 0
                    CalcCnt = CalcCnt + 1;
                    if CalcCnt > 1
                        disp(ptlin); disp(ptskip)
                        if fidO > 0
                            fprintf(fidO, ptlin); fprintf(fidO, ptskip);
                        end
                    end
                    pt2 = sprintf('X%0.f  %s      X%.0f  %s \n', ...
                           iii1, xxlbl, iii2, yylbl);
                           
                    disp(pt3); disp(pt2); disp(ptskip)
                    if fidO > 0 
                        fprintf(fidO, pt3); fprintf(fidO, pt2);
                        fprintf(fidO, ptskip);
                    end
                  if dothis == 3
                      RtrnCode = Vector_Corr_PMCorr(LV, LV1, alfa, fidO);
                  else
                    RtrnCode = Vector_Corr_LinCirc1(AV1, LV, alfa, fidO, ...
                                                    NTrials);
                    RtrnCode = Vector_Corr_CAssoc1(AV1, LV, alfa, fidO, ...
                                                   NTrials);
                    % following converts linear to circular, then processes
                    % Fisher, p. 161
                    AzL = 2*atan(LV);
                    RtrnCode = Vector_Corr_CircTlin1(AV1, AzL, alfa, ...
                                                    fidO, 1, NTrials);
                    if Plots > 0
                      [RtrnCode, figCntX] = Vector_Corr_Plots(2, ...
                          AV1, LV, xlbl, ylbl, DataTtl, figCnt);
                      figCnt = figCntX;
                    end
                  end  
                end       
            end
        end
    end
end

 
