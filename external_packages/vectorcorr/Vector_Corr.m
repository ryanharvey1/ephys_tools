%Vector_Corr.m

% This script is described in the Computers & Geosciences paper:
%   MATLAB functions to analyze directional (azimuthal) data.
%   II. Correlation
%   by Thomas A. Jones
%   2006, v. 32

% Copyright C 2004, 2005  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Calculates linear-circular and circular-circular correlation analyses
% for a dataset containing a mixture of azimuth and ordinary variables.  
% The vectors (azimuths) input are in degrees, defined over (0,360)
% [equivalent to (0, 2PI) radians].

% User controls primarily are specified via a GUI.

% Vector and linear data may be read from an external ascii file, 
% or in a 2D array already in Matlab memory.  For both, Rows represent 
% various data values (observations), and two or more of the columns  
% contain the Azimuth and/or linear data.   

% Calculations include:
%   Cross plots of the azimuth-azimuth pairs or of the linear-azimuth pairs
%   Azimuth-azimuth and linear-azimuth correlations, and tests of 
%     significance where sample size is great enough
%   Simple statistics (azimuthal or linear) on each variable
%   No linear-linear correlations are calculated

% Commonly used variables in the functions include:
%   AzData: input array containing Azimuth and independent variable(s)
%   IndepX: list of which columns in AzData are used as indep. variables
%   IndepName: identifiers of independent variables
%   IndepAz: indicates if variable is linear(=0) or azimuth (=1)
%   Xarray: array of variables selected to use 
%   NTot: N, sample size
%   DataTtl: job title
%   alfa: significance level for tests
%   fidO: output file for summary of job (>0 implies write file)
%   Plots: are plots to be generated (=1) 
%   RtrnCode: error return code
%    0 = OK
%    1 = Input data set could not be opened or not found
%    2 = Incorrect numbers of columns specified for input dataset,
%        or it contains characters/strings, or incomplete record.
%    3 = Could not open external output file for calculations
%    7 = User didn't specify name of external output file
%    8 = Job cancelled while specifying job/data controls

% Functions and scripts called from this module:
%   Vector_DataTypeG.m
%   Vector_Corr_CntlsG_input.m
%   Read_Reals_File.m
%   Vector_Corr_Calcs.m

% set up for job and printing

pthdr=sprintf('CORRELATION ANALYSIS FOR VECTORIAL DATA\n');
ptskip=sprintf('\n');
pteq=sprintf('==================================================\n');

disp(ptskip); disp(pteq); disp(ptskip); disp(pthdr); disp(ptskip); 

RtrnCode = 0;
fclose all;

% determine if using data array or input external file
% get job title or description

datatype = 0;
DataTtl = '';
Vector_DataTypeG;

disp(DataTtl);

% Get name of external file to analyze

if datatype == 1
   disp(' '); disp(' ')
   disp('Specify external text file containing data')

   [fileIn,pathIn] = uigetfile('*.txt','Files of type txt',10,10);
   if ~isstr(fileIn);
       disp('No external file specified')
       RtrnCode=7;
       ptstop=sprintf('JOB TERMINATING: ReturnCode = %.0f \n',RtrnCode);
       disp(ptskip); disp(ptstop); disp(pteq); disp(ptskip)
       return
   end
end

% Input general controls on data set to be used  
    
Vector_Corr_CntlsG_input;
if RtrnCode > 0
    ptstop=sprintf('JOB TERMINATING: ReturnCode = %.0f \n',RtrnCode);
    disp(ptskip); disp(ptstop); disp(pteq); disp(ptskip)
    return
end
NAzimVars = cntAz;

% Input data array or external file

if datatype == 0
    AzData = eval(Data_Array);
    ptdsnIn=['Name of input data array: ', Data_Array]; 
    disp(ptdsnIn);   disp(ptskip)
else
   dsnIn = [pathIn, fileIn];
   [RtrnCode, NRows, AzData] = Read_Reals_File(dsnIn, NSkipHdr, NCols);
   if RtrnCode > 0
      ptstop=sprintf('JOB TERMINATING: ReturnCode = %.0f \n',RtrnCode);
      disp(ptskip); disp(ptstop); disp(pteq); disp(ptskip)
      return
   end
   ptdsnIn=['Name of input data file: ', fileIn];
   disp(ptdsnIn);   disp(ptskip)
end

% open file to write summaries

fidO = -99; 
if WrOut == 0
    dsnOut = [pathOut, fileOut];
    fidO = fopen(dsnOut, 'wt');
    if fidO < 0 
        RtrnCode = 3;
        disp('ERROR *** Output file could not be opened')
        ptstop=sprintf('JOB TERMINATING: ReturnCode = %.0f \n',RtrnCode);
        disp(ptskip); disp(ptstop); disp(pteq); disp(ptskip)
        return
    end
    disp(ptskip)
    ptdsnOut=['File opened for output of job summary: ', fileOut];
    
    % write job headers info to file

    fprintf(fidO, ptskip);   fprintf(fidO, pthdr);
    fprintf(fidO, ptskip);
    
    fprintf(fidO, [DataTtl,'\n']); fprintf(fidO, ptskip);
    fprintf(fidO, ptdsnIn);        fprintf(fidO, ptskip); 
    fprintf(fidO, ptskip);
end

% summarize information about data set

pt1=sprintf('Data array contains %.0f Rows and %.0f Columns\n', ...
             NRows, NCols);
if DegRad == 0
    pt1a = sprintf('Azimuths input as Degrees\n');
else
    pt1a = sprintf('Azimuths input as Radians\n');
end
pt2 = sprintf(['Sample size: N = %.0f\n' ...
               'Significance level input: alfa = %.3f\n'], NRows, alfa);
pt2a = sprintf('Number of trials used for resampling: %.0f\n', NTrials);       
pt3=sprintf('Variables     col.       Azim?\n             in data\n');         
disp(pt1); disp(pt1a); disp(pt2); 
if NRows < 36
    disp(pt2a)
end
disp(ptskip); disp(pt3);
for iii = 1:MaxIndep
    if IndepX(iii) > 0
        if IndepAz(iii) > 0
            Az = 'Yes';
        else
            Az = 'No ';
        end
        pt4=sprintf('    X%.0f          %.0f         %s     %s', ...
            iii, IndepX(iii), Az, IndepName(iii, :));
    else
        pt4=sprintf('    X%.0f:  Not specified', iii);
    end
    disp(pt4)
end

if fidO > 0
    fprintf(fidO, pt1);      fprintf(fidO, pt1a);   fprintf(fidO, pt2); 
    if NRows < 36
       fprintf(fidO, pt2a);
    end
    fprintf(fidO, ptskip);   fprintf(fidO, pt3);
    for iii = 1:MaxIndep
     if IndepX(iii) > 0
        if IndepAz(iii) > 0
            Az = 'Yes';
        else
            Az = 'No ';
        end
        pt4=sprintf('    X%.0f          %.0f         %s     %s\n', ...
            iii, IndepX(iii), Az, IndepName(iii, :));
     else
        pt4=sprintf('    X%.0f:  Not specified\n', iii); 
     end
     fprintf(fidO, pt4);
    end
    
    fprintf(fidO, ptskip);fprintf(fidO, ptdsnOut);
    disp(ptskip); disp(ptdsnOut); disp(ptskip); disp(pteq)
end
    
% Process the data - prepare and do correlation analyses

RtrnCode = Vector_Corr_Calcs(AzData, DegRad, NIndepX, IndepX, ... 
    IndepName, IndepAz, alfa, fidO, Plots, DataTtl, NTrials);

disp(' ');
if RtrnCode == 0 
    disp('JOB COMPLETED');
end
disp('==================================================')
disp(' ')
if fidO > 0
    fclose(fidO);
end

