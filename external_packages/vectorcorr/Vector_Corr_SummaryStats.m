function RtrnCode = Vector_Corr_SummaryStats(DataArray, fidO, ...
             IndepX, IndepName, IndepAz)

% Vector_Corr_SummaryStats.m

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Calculate summary statistics for the variables
%   For azimuth variables, vector mean and concentration
%   For linear variables, simple mean, variance, etc

% Input variables:
%   DataArray: array of variables selected for analysis (degrees)
%   fidO: indicates file to output summary stats to (if > 0)
%   IndepX: list of variables to process 
%   IndepName: names of variables
%   IndepAz: which variables are azimuth (=1)
% Output variables:
%   RtrnCode:  Error return code 

% Functions called by this module:
%   VectMean_arctan.m
%   CalcKappa.m


RtrnCode = 0;
[NTot, NC] = size(DataArray);
ycos = zeros(1,NTot);  ysin = zeros(1,NTot);

ptskip = sprintf('\n');
pteq = sprintf('====================================================\n');
ptlin = sprintf('----------------------------------------------------\n');
disp(ptskip); disp(pteq); 
if fidO > 0
    fprintf(fidO, ptskip); fprintf(fidO, pteq); 
end

jjj = 0;
for iii = 1:12
    if IndepX(iii) > 0
        jjj = jjj + 1;
        
        % put header for this variable
        pt1 = sprintf('X%.0f    col. %.0f    %s\n', ...
              jjj, IndepX(iii), IndepName(iii,:));
        disp(ptskip)
        if jjj > 1; disp(ptlin); end
        disp(ptskip); disp(pt1)
        if fidO > 0
            fprintf(fidO, ptskip);
            if jjj == 1; fprintf(fidO, ptlin); end
            fprintf(fidO, ptskip); fprintf(fidO, pt1);
        end
        
        if IndepAz(iii) > 0
        
            % Azimuths
            ycos = cos(DataArray(:, jjj)/57.3);
            ysin = sin(DataArray(:, jjj)/57.3);          
            Uval = sum(ycos);
            Vval = sum(ysin);
            R2val = Uval^2 + Vval^2;
            Rval = sqrt(R2val); 
            Rbar = Rval/NTot;
            ThetaHat = VectMean_arctan(Vval, Uval);
            KappaHat = CalcKappa(Rbar, 0);
            if NTot < 20
              if KappaHat < 2
                 KappaHatCorr = max(KappaHat - 2/(NTot*KappaHat), 0);
              else
                 KappaHatCorr = KappaHat*(NTot - 1)^3 / (NTot + NTot^3);
              end
            else
                 KappaHatCorr = -99;
            end

            pt21 = sprintf(['     R-square = %9.5g \n', ...
                    '     R-bar = %9.5g\n', ...
                    '     VectorMean (deg) = %6.1f\n', ...
                    '     Kappa-hat = %9.5g\n'], ...
                    R2val, Rbar, ThetaHat, KappaHat); 
            disp(pt21)
            if KappaHatCorr > -5
             pt6=sprintf(['     Kappa-hat (corrected for small-sample ',...
                          'bias) = %.5g\n', ...
                          '     (Ref.: Fisher, 1993, p. 88 (4.41))\n'],...
                           KappaHatCorr);             
             disp(pt6) 
            end
            
            if fidO > 0
               fprintf(fidO, ptskip);
               fprintf(fidO, pt21); 
               if KappaHatCorr > -5 
                  fprintf(fidO, [pt6,'\n']); 
               end
            end
              
        else
            
            % linear variables
           Xmn = mean(DataArray(:,jjj));
           Xvar = var(DataArray(:,jjj));  Xsd = sqrt(Xvar);
           pt31 = sprintf(['     Mean = %9.5g\n     Variance = %9.5g\n',...
               '     Standard deviation = %9.5g \n'], Xmn, Xvar, Xsd);
           disp(pt31);
           if fidO > 0
               fprintf(fidO, pt31);
           end
       end
   end
end
        

