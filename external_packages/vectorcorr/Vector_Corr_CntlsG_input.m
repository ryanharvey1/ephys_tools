%Vector_Corr_CntlsG_input.m

% Copyright C 2004  Thomas A. Jones.  All rights reserved.
% Everyone is granted permission to copy, modify, and redistribute
% this program, but under the condition that this copyright notice
% remain intact and the Computers & Geosciences source be referenced. 
% This software has been tested extensively, but is distributed 
% WITHOUT ANY WARRANTY.  The author does not accept responsibility
% for any possible consequences of using it.

% Driver for GUI that obtains the controls that specify/define
% the dataset to be used (either external file or existing data array),
% and which variables in dataset for linear-circular and 
% circular-circular correlation analysis. 
% Also gets controls for calculations.
% This drives the obtaining of values and checking for errors.

% Variables to be obtained and defaults:
%   NSkipHdr=0: Number of headers on data file to be skipped
%   NCols=0: Number of columns in data file
%   X1,X2,X3,...,X12: Which columns in data file/set contains X variables
%   XName1, ..., XName12: Alpha identifiers of X variables
%   Az1 ... Az12: Which columns are azimuths? (=1)
%   DegRad: Input degrees (0) or radians(1)
%   NTrials: Number of resampling trials (100 - 10000) default=2500
%   alfa: significance level for tests 
%     should be 0.10, 0.05, or 0.01 for CAssociation test
%     should be 0.10, 0.05, 0.025, 0.01, 0.005 for CircularRank test
%   Plots = 0: Are crossplots to be made? (yes=1)
%   WrOut = 1: Write summary of job to external file (1=No, 0=Yes)
%   Data_Array='': Name of array that contains data 
%      default = 'AzData' if currently in memory
%   RtrnCode = 0: Error return

% Script called from this module:
%   Vector_Corr_CntlsG.m

% initialize variables to be obtained from GUI

NSkipHdr=0;
NCols=0;
DegRad = 0;
NTrials = 2500;
X1=0; X2=0; X3=0; X4=0; X5=0; X6=0; X7=0; X8=0; X9=0; X10=0;
X11=0; X12=0;
XName1=''; XName2=''; XName3=''; XName4=''; XName5='';
XName6=''; XName7=''; XName8=''; XName9=''; XName10='';
XName11=' '; XName12=' '; 
Az1=0; Az2=0;Az3=0;Az4=0;Az5=0;Az6=0;Az7=0;Az8=0;Az9=0;Az10=0;
Az11=0;Az12=0;
Plots = 0;
alfa = 0.05;
WrOut = 1;
Data_Array=''; 
if datatype ~= 1
    slist = who;
    for iiii = 1:length(slist)
        if strcmp(slist(iiii), 'AzData');
            Data_Array = 'AzData';
            break
        end
    end
end
RtrnCode = 0;
IndepX = zeros(1,12); IndepAz = zeros(1,12); 

% Loop over the GUI figure.  For each iteration,
% obtain values, check for errors, notify, re-call GUI

done=0;
while done==0
    
    % Get controls  
    
    Vector_Corr_CntlsG;
    
    disp(' '); disp('-----------------------------------------------')
    disp(' ')
    
    if RtrnCode == 8
        disp('Job cancelled by user')
        done = 1;
        return
    end
    
    %Check for errors/omissions
    
    done = 1;

  if datatype == 1
        
       % Header skips
    if NSkipHdr < 0 
        done = 0;
        prt1=['ERROR *** Invalid number of headers to skip: ', ...
                num2str(NSkipHdr)];
        disp(prt1)
    end
    
       % Number columns in file
    if NCols < 2
        done = 0;
        prt1=['ERROR *** Invalid number of columns in datafile: ', ...
                num2str(NCols)];
        disp(prt1)
    end
    
  else
    
    if ~isstr(Data_Array); 
        done = 0;
        disp('ERROR *** No name specified for data array to process')
    else 
        NRows = 0; NCols = 0;
        slist = who;
        for iiii=1:length(slist);
            if strcmp(slist(iiii), Data_Array);
               aaa = eval(Data_Array);
               [NRows,NCols] = size(aaa);
               break
            end
        end
            if NRows == 0
               done = 0;
               prt1=['ERROR *** Data-array with specified name does',...
                       ' not exist: ',Data_Array];
               disp(prt1)
            end
    end
    
end
    
       % alfa    
     if alfa > 0.25 
        done = 0;
        prt1=['ERROR *** Invalid significance level (alfa): ', ...
                num2str(alfa)];
        disp(prt1)
        disp( 'Allowable value is in range (0.001, 0.25)') 
     end

            % NTrials
     if NTrials < 100 | NTrials > 10000 
        done = 0;
        prt1=['ERROR *** Invalid number of resampling trials: ', ...
                num2str(NTrials)];
        disp(prt1)
        disp( 'Allowable value is in range (100, 10000)') 
     end
    
    NIndepX = 0; MaxIndep = 0;
    IndepX(1) = X1; IndepX(2) = X2; IndepX(3) = X3;
    IndepX(4) = X4; IndepX(5) = X5;
    IndepX(6) = X6; IndepX(7) = X7; IndepX(8) = X8;
    IndepX(9) = X9; IndepX(10) = X10; IndepX(11) = X11; IndepX(12) = X12;
    IndepName = str2mat(XName1, XName2, XName3, XName4, XName5, ...
                XName6, XName7, XName8, XName9, XName10, XName11, XName12);
    IndepAz(1) = Az1; IndepAz(2) = Az2; IndepAz(3) = Az3; 
    IndepAz(4) = Az4; IndepAz(5) = Az5; IndepAz(6) = Az6;
    IndepAz(7) = Az7; IndepAz(8) = Az8; IndepAz(9) = Az9;
    IndepAz(10) = Az10; IndepAz(11) = Az11; IndepAz(12) = Az12;
    
    % check inputs of variables
    
    for iii = 1:12  
        if IndepX(iii) ~= 0
            NIndepX = NIndepX + 1;
            MaxIndep = iii;
            
          % check validity of variables
           if IndepX(iii) > NCols
             done = 0;
             prt1=['ERROR *** Invalid column number for X',num2str(iii),...
                   ' data: ', num2str(IndepX(iii)), ' (', ...
                   IndepName(iii,:), ')'];
             disp(prt1);
             prt1=['          Col. specified > max columns in data: ',...
                   num2str(NCols)];
             disp(prt1);
          end
          if IndepX(iii) < 0
             done = 0;
             prt1=['ERROR *** Invalid column number for X',num2str(iii),...
                   ' data: ', num2str(IndepX(iii)), ' (', ...
                   IndepName(iii,:), ')'];
             disp(prt1);
          end
      end
  end
    
     % check if any azimuths specified
     
     cntAz = 0;
     for iii = 1:12  
        if IndepX(iii) > 0
            if IndepAz(iii) > 0 
                cntAz = cntAz + 1;
            end
        end
    end
    if cntAz == 0
             done = 0;
             prt1='ERROR *** No Azimuth variables were specified';  
             disp(prt1);
    end

  % check if named no variables or if any numbers duplicated
  if NIndepX == 0
      done = 0;
      disp('ERROR *** No variables were specified');
  else
      for iii = 1:MaxIndep - 1
          for jjj = iii+1:MaxIndep
              if (IndepX(iii) == IndepX(jjj)) & (IndepX(iii) > 0) 
                  done = 0;
                  prt1=['ERROR *** Variables X%0.f and', ...
                        ' X%.0f specified with same column number'];
                  prt2=sprintf(prt1, iii, jjj);
                  disp(prt2);
              end
          end
      end
  end
  
  if done == 0
     disp(' '); 
     prt1=['Number of variables specified: ',num2str(NIndepX)];
     disp(prt1);
      for iii = 1:MaxIndep
          prt1=['    X', num2str(iii), '      ',num2str(IndepX(iii)), ...
               '       ', IndepName(iii,:)];
          disp(prt1);
      end
      disp(' ');
  end
  
end

% specify name of output file

if WrOut == 0
    
    [fileOut, pathOut] = uiputfile('*.txt','Specify output file',10, 10);
    
    if ~isstr(fileOut);
        disp('No output file specified by user - Job cancelled')
        RtrnCode = 7;
        return
    end
    
end
