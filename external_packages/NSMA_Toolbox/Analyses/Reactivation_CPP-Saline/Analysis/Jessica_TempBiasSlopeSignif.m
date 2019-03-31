% script to calculate significance of the slope differences in CPP and Saline TempBias Regression analysis
%
% PL April 11, 2004

% run  mat files for CPP and Saline regression results

runNrStr = 'Run000_'
outDir = 'T:\CPP1\JessicasData\Analysis\out'
groupName = {'Saline', 'CPP'};

reRun = 1;   % set to 1 if you are re-running the script for the second part and have already all variables from 1st part in workspace

diary off;                  % in case it is on from a bombed previous run
timestamp = strrep(datestr(datenum(now),0),':','-');
diary([outDir, filesep, 'CombinedRegressAnalysisDiary-' , 'run_',timestamp,'.txt']);


plotRats = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1st part: compute all regression slopes and regress stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over saline and CPP
if ~reRun

clear res;
    
for CPP = 0:1
        
	if ~CPP
        run_label = 'Saline';
        indir = [outDir '\TempBias_Saline_15min'];
        ratlist = {'7508_200ms';'7651_200ms';'7653_200ms';'7658_200ms'};
	else    
        run_label = 'CPP';
        indir = [outDir '\TempBias_CPP_15min'];
        ratlist = {'7508_200ms';'7651_200ms';'7653_200ms';'7658_200ms'};
	end
	
	
	outputDir = indir
	outPrefix = [runNrStr 'RegAnalysis_'];
	
	
	%for diary output
	run_label
	CPP
	outputDir
	ratlist
	
	
	%initialize
	clear S*; clear B*; clear sys*; clear A*; clear X;
	S1_bias_pooled = [];
	B_bias_pooled = [];
	S2_bias_pooled = [];
	
	%for ii = 1:length(bias_list)
	for jj = 1:length(ratlist)  
        pushdir(indir);
        pwd;
        
        
        filename = [runNrStr 'Bias_results_' ratlist{jj} '.mat'];
        load(filename);
        fprintf('\n\n=========================== %s ============================\n', filename);
        
        X = [ones(size(ABM,1),1) ABM(:,2)];     % ABM is the "all bias matrix which has three columns representing all the bias values for  % S1 - B - S2 respectively
        B_bias_pooled = [B_bias_pooled; ABM(:,2)];
        S1_bias_pooled = [S1_bias_pooled; ABM(:,1)];
        if size(ABM,1) > 1
            [S1Bstats.b,S1Bstats.bint,S1Bstats.r,S1Bstats.rint,S1Bstats.stats] = regress(ABM(:,1),X);  % S1-B 
        else
            S1Bstats.b = NaN;
            S1Bstats.bint = NaN;
            S1Bstats.stats = NaN;
        end
            
        disp('S1Bstats.b ='); S1Bstats.b
        disp('S1Bstats.bint ='); S1Bstats.bint
        disp('S1Bstats.stats ='); S1Bstats.stats
        
        res{CPP+1}.x{jj} = ABM(:,2);
        res{CPP+1}.S1Bstats{jj}.b = S1Bstats.b;
        res{CPP+1}.S1Bstats{jj}.r = S1Bstats.r;
        res{CPP+1}.S1Bstats{jj}.bint = S1Bstats.bint;
        res{CPP+1}.S1Bstats{jj}.stats = S1Bstats.stats;
        
        
        if plotRats & (size(ABM,1) > 1)
            fh = figure;
			scrPos = get(0,'ScreenSize');
			set(fh,'Units','pixels');
			set(fh,'Position',[0,0,(scrPos(4)-20)/2,scrPos(4)-20]);
            orient tall;
            subplot(2,1,1);
            plot(ABM(:,2),ABM(:,1),'.') 
            axis equal
			hold on
			reg = S1Bstats;
			x = [-1.5,1.5];
			y = reg.b(2)*x + reg.b(1);
			lh = line(x,y);
			set(lh,'LineStyle','-','Color','red');
			textstr = sprintf('%6.3f*x %+6.3f',reg.b(2),reg.b(1));
			text(-1,0.9,textstr);
			hold off
            grid on
            [NDots,Ncols] = size(ABM)
            title([runNrStr ' ' ratlist{jj} ' ' run_label ': B vs S1    (' num2str(NDots) ' data points)'])
            xlabel('Behavior')
            ylabel('Sleep 1')
            axis([-1,1,-1,1])
            axis equal
        end
        
        if size(ABM,1) > 1
            [S2Bstats.b,S2Bstats.bint,S2Bstats.r,S2Bstats.rint,S2Bstats.stats] = regress(ABM(:,3),X);  % S2-B
         else
            S2Bstats.b = NaN;
            S2Bstats.bint = NaN;
            S2Bstats.stats = NaN;
        end
        S2_bias_pooled = [S2_bias_pooled; ABM(:,3)];
        disp('S2Bstats.b ='); S2Bstats.b
        disp('S2Bstats.bint ='); S2Bstats.bint
        disp('S2Bstats.stats ='); S2Bstats.stats
        
        res{CPP+1}.S2Bstats{jj}.b = S2Bstats.b;
        res{CPP+1}.S2Bstats{jj}.r = S2Bstats.r;
        res{CPP+1}.S2Bstats{jj}.bint = S2Bstats.bint;
        res{CPP+1}.S2Bstats{jj}.stats = S2Bstats.stats;
        
        if plotRats & (size(ABM,1) > 1)
            subplot(2,1,2);
            plot(ABM(:,2),ABM(:,3),'.') 
			hold on
			reg = S2Bstats;
			x = [-1.5,1.5];
			y = reg.b(2)*x + reg.b(1);
			lh = line(x,y);
			set(lh,'LineStyle','-','Color','red');
			textstr = sprintf('%6.3f*x %+6.3f',reg.b(2),reg.b(1));
			text(-1,0.9,textstr);
			hold off
            grid on
            title([runNrStr ' ' ratlist{jj} ' ' run_label ': B vs S2    (' num2str(NDots) ' data points)'])
            xlabel('Behavior')
            ylabel('Sleep 2')
            axis([-1,1,-1,1])
            axis equal
        end
        
        
        fprintf([' S2B_allbiaspart = %f; S2B_allbiaspart_rev = %f; \n' ...
                '  S1B_allbiascorr = %f; S2B_allbiascorr = %f\n\n'],S2B_allbiaspart,S2B_allbiaspart_rev,S1B_allbiascorr,S2B_allbiascorr);
        
        s1file = [outputDir, filesep, outPrefix, run_label, '_', ratlist{jj}, '_S1Bstats'];
        s2file = [outputDir, filesep, outPrefix, run_label, '_', ratlist{jj}, '_S2Bstats'];    
        save (s1file, 'S1Bstats');
        save (s2file, 'S2Bstats');
        if plotRats
            saveas(gcf,[outputDir, filesep, outPrefix, run_label, '_', ratlist{jj}],'fig');
        end
            
        %clear S*; clear B*; clear sys*; clear A*; clear X;
        
        popdir;
        fprintf('======================================================================\n');
        
	end %for jj


    
	fprintf('======================================================================\n');
	fprintf('======================================================================\n');
	fprintf('============ All Rats All Session pooled together :===================\n\n');
	[NDots,Ncols] = size(S1_bias_pooled);
	NDots 

    res{CPP+1}.B_bias_pooled = B_bias_pooled;
    res{CPP+1}.S1_bias_pooled = S1_bias_pooled;
    res{CPP+1}.S2_bias_pooled = S2_bias_pooled;
    
	X = [ones(size(B_bias_pooled,1),1) B_bias_pooled];
	[S1Bstats.b,S1Bstats.bint,S1Bstats.r,S1Bstats.rint,S1Bstats.stats] = regress(S1_bias_pooled,X);  % S1-B 
	disp('S1Bstats.b ='); S1Bstats.b
	disp('S1Bstats.bint ='); S1Bstats.bint
	disp('S1Bstats.stats ='); S1Bstats.stats
    
    res{CPP+1}.S1BstatsPooled.b = S1Bstats.b;
    res{CPP+1}.S1BstatsPooled.r = S1Bstats.r;
    res{CPP+1}.S1BstatsPooled.bint = S1Bstats.bint;
    res{CPP+1}.S1BstatsPooled.stats = S1Bstats.stats;

	fh = figure;
	scrPos = get(0,'ScreenSize');
	set(fh,'Units','pixels');
	set(fh,'Position',[0,0,(scrPos(4)-20)/2,scrPos(4)-20]);
	orient tall;
	subplot(2,1,1);
	plot(B_bias_pooled,S1_bias_pooled,'.') 
	hold on
	reg = S1Bstats;
	x = [-1.5,1.5];
	y = reg.b(2)*x + reg.b(1);
	lh = line(x,y);
	set(lh,'LineStyle','-','Color','red');
	textstr = sprintf('%6.3f*x %+6.3f',reg.b(2),reg.b(1));
	text(-1,0.9,textstr);
	hold off
	grid on
	title([runNrStr ' ' run_label ': ' num2str(length(ratlist)) 'Rats combined: B vs S1    (' num2str(NDots) ' data points)'])
	xlabel('Behavior')
	ylabel('Sleep 1')
	axis([-1,1,-1,1])
	axis equal
	
	
	[S2Bstats.b,S2Bstats.bint,S2Bstats.r,S2Bstats.rint,S2Bstats.stats] = regress(S2_bias_pooled,X);  % S2-B 
	disp('S2Bstats.b ='); S2Bstats.b
	disp('S2Bstats.bint ='); S2Bstats.bint
	disp('S2Bstats.stats ='); S2Bstats.stats

    res{CPP+1}.S2BstatsPooled.b = S2Bstats.b;
    res{CPP+1}.S2BstatsPooled.r = S2Bstats.r;
    res{CPP+1}.S2BstatsPooled.bint = S2Bstats.bint;
    res{CPP+1}.S2BstatsPooled.stats = S2Bstats.stats;

	subplot(2,1,2);
	plot(B_bias_pooled,S2_bias_pooled,'.') 
	hold on
	reg = S2Bstats;
	x = [-1.5,1.5];
	y = reg.b(2)*x + reg.b(1);
	lh = line(x,y);
	set(lh,'LineStyle','-','Color','red');
	textstr = sprintf('%6.3f*x %+6.3f',reg.b(2),reg.b(1));
	text(-1,0.9,textstr);
	hold off
	grid on
	title([runNrStr ' ' run_label ': ' num2str(length(ratlist)) ' Rats combined: B vs S2    (' num2str(NDots) ' data points)'])
	xlabel('Behavior')
	ylabel('Sleep 2')
	axis([-1,1,-1,1])
	axis equal
	
	saveas(gcf,[outputDir, filesep, outPrefix, run_label,'_AllRatsCombined'],'fig');
	
	allfile = [outputDir, filesep, outPrefix, run_label, '_', 'AllRatsCombined_Allstats'];
	save(allfile,'*_pooled','S1Bstats','S2Bstats');
	fprintf('================= Done ================================================\n');
	
end % for CPP

end % if rerun

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2nd part: significance tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Run sigifcance test for slope differences for each rat and all rats pooled 

filen = [outputDir, filesep, outPrefix, run_label, '_', 'AllstatsForSignificanceTest'];
save(filen,'res');

% open csv output file
csvfile = [outDir, filesep, 'CombinedRegressAnalysis-AllRatsCombined_T-Tests-' , 'run_', timestamp,'.csv'];
fpCSV = fopen(csvfile ,'wt');
% write CSV headers
fprintf(fpCSV,['Rat,','TestType,', 'Slope 1,', 'Slope 2,','sig1,','sig2,','#DOF 1,', '#DOF2,','alpha,','tail,','t,','p,','h,','s1-s2,','ci_lo,','ci_hi,' ,'\n\n']);

for jj = 1:length(ratlist)
    % extract s1, s2, sig1, sig2 and n1, n2  For S1B slopes in saline and CPP
    s1 = res{1}.S1Bstats{jj}.b(2);
    x1 = res{1}.x{jj};
    n1 = length(x1);
    r1 = res{1}.S1Bstats{jj}.r;
    sig1 = sqrt( sum(r1.^2)/((n1-2)*sum((x1-mean(x1)).^2)) );  

    s2 = res{2}.S1Bstats{jj}.b(2);
    x2 = res{2}.x{jj};
    n2 = length(x2);
    r2 = res{2}.S1Bstats{jj}.r;
    sig2 = sqrt( sum(r2.^2)/((n2-2)*sum((x2-mean(x2)).^2)) );  
    
    
    % extract s3, s4, sig3, sig4 and n3, n4  For S2B slopes in saline and CPP
    s3 = res{1}.S2Bstats{jj}.b(2);
    x3 = res{1}.x{jj};
    n3 = length(x1);
    r3 = res{1}.S2Bstats{jj}.r;
    sig3 = sqrt( sum(r3.^2)/((n3-2)*sum((x3-mean(x3)).^2)) );  

    s4 = res{2}.S2Bstats{jj}.b(2);
    x4 = res{2}.x{jj};
    n4 = length(x2);
    r4 = res{2}.S2Bstats{jj}.r;
    sig4 = sqrt( sum(r4.^2)/((n4-2)*sum((x4-mean(x4)).^2)) );  
    
    
    % compute significances between S1B and S2B slopes
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    alpha = 0.05
    fprintf(fpCSV,'\n\n');
    disp( sprintf('\n\n =======   rat : %s  S1B(Saline) and S2B(Saline) significance test:', ratlist{jj}) );
    tail = -1
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s3,sig3,n3, alpha, tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S1B(Saline) and S2B(Saline)',s1,s3,sig1,sig3,n1,n3,alpha,tail,t,p,h,s1-s3,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);
    tail = +1
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s3,sig3,n3, alpha, tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S1B(Saline) and S2B(Saline)',s1,s3,sig1,sig3,n1,n3,alpha,tail,t,p,h,s1-s3,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);
    tail = 0
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s3,sig3,n3, alpha, tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S1B(Saline) and S2B(Saline)',s1,s3,sig1,sig3,n1,n3,alpha,tail,t,p,h,s1-s3,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);
    
    fprintf(fpCSV,'\n');
    disp( sprintf('\n\n =======   rat : %s  S1B(CPP) and S2B(CPP) significance test:', ratlist{jj}) );
    tail = -1
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s2,sig2,n2,s4,sig4,n4, alpha, tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S1B(CPP) and S2B(CPP)',s2,s4,sig2,sig4,n2,n4,alpha,tail,t,p,h,s2-s4,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);
    tail = -0
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s2,sig2,n2,s4,sig4,n4, alpha, tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S1B(CPP) and S2B(CPP)',s2,s4,sig2,sig4,n2,n4,alpha,tail,t,p,h,s2-s4,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);
    tail = -0
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s2,sig2,n2,s4,sig4,n4, alpha, tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S1B(CPP) and S2B(CPP)',s2,s4,sig2,sig4,n2,n4,alpha,tail,t,p,h,s2-s4,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);

    fprintf(fpCSV,'\n');
    disp( sprintf('\n\n =======   rat : %s  S1B(saline) and S1B(CPP) significance test:', ratlist{jj}) );
    tail = -1
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s2,sig2,n2, alpha,tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S1B(saline) and S1B(CPP)',s1,s2,sig1,sig2,n1,n2,alpha,tail,t,p,h,s1-s2,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);
    tail = +1
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s2,sig2,n2, alpha,tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S1B(saline) and S1B(CPP)',s1,s2,sig1,sig2,n1,n2,alpha,tail,t,p,h,s1-s2,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);
    tail = 0
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s2,sig2,n2, alpha,tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S1B(saline) and S1B(CPP)',s1,s2,sig1,sig2,n1,n2,alpha,tail,t,p,h,s1-s2,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);

    
    fprintf(fpCSV,'\n');
    disp( sprintf('\n\n =======   rat : %s  S2B(saline) and S2B(CPP) significance test:', ratlist{jj}) );
    tail = -1
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s3,sig3,n3,s4,sig4,n4, alpha,tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S2B(saline) and S2B(CPP)',s3,s4,sig3,sig4,n3,n4,alpha,tail,t,p,h,s3-s4,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);
    tail = +1
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s3,sig3,n3,s4,sig4,n4, alpha,tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S2B(saline) and S2B(CPP)',s3,s4,sig3,sig4,n3,n4,alpha,tail,t,p,h,s3-s4,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);
    tail = 0
    [t, p, ci, ndof, h] = SlopeDiff_TTest(s3,sig3,n3,s4,sig4,n4, alpha,tail)
    cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
        ratlist{jj},'S2B(saline) and S2B(CPP)',s3,s4,sig3,sig4,n3,n4,alpha,tail,t,p,h,s3-s4,ci(1),ci(2));
    fprintf(fpCSV,cvsstr);
    
end % for jj

%%%%%%%%%%%%%%   all rats combined

% extract s1, s2, sig1, sig2 and n1, n2  For S1B slopes in saline and CPP
s1 = res{1}.S1BstatsPooled.b(2);
x1 = res{1}.B_bias_pooled;
n1 = length(x1);
r1 = res{1}.S1BstatsPooled.r;
sig1 = sqrt( sum(r1.^2)/((n1-2)*sum((x1-mean(x1)).^2)) );  

s2 = res{2}.S1BstatsPooled.b(2);
x2 = res{2}.B_bias_pooled;
n2 = length(x2);
r2 = res{2}.S1BstatsPooled.r;
sig2 = sqrt( sum(r2.^2)/((n2-2)*sum((x2-mean(x2)).^2)) );  


% extract s3, s4, sig3, sig4 and n3, n4  For S2B slopes in saline and CPP
s3 = res{1}.S2BstatsPooled.b(2);
x3 = res{1}.B_bias_pooled;
n3 = length(x1);
r3 = res{1}.S2BstatsPooled.r;
sig3 = sqrt( sum(r3.^2)/((n3-2)*sum((x3-mean(x3)).^2)) );  

s4 = res{2}.S2BstatsPooled.b(2);
x4 = res{2}.B_bias_pooled;
n4 = length(x2);
r4 = res{2}.S2BstatsPooled.r;
sig4 = sqrt( sum(r4.^2)/((n4-2)*sum((x4-mean(x4)).^2)) );  


% compute significances between S1B and S2B slopes
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALL RATS COMBINED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
alpha = 0.05
fprintf(fpCSV,'\n\n');
disp( sprintf('\n\n =======   All rats pooled:  S1B(Saline) and S2B(Saline) significance test:') );
tail = -1
[t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s3,sig3,n3, alpha, tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S1B(Saline) and S2B(Saline)',s1,s3,sig1,sig3,n1,n3,alpha,tail,t,p,h,s1-s3,ci(1),ci(2));
fprintf(fpCSV,cvsstr);
tail = +1
[t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s3,sig3,n3, alpha, tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S1B(Saline) and S2B(Saline)',s1,s3,sig1,sig3,n1,n3,alpha,tail,t,p,h,s1-s3,ci(1),ci(2));
fprintf(fpCSV,cvsstr);
tail = 0
[t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s3,sig3,n3, alpha, tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S1B(Saline) and S2B(Saline)',s1,s3,sig1,sig3,n1,n3,alpha,tail,t,p,h,s1-s3,ci(1),ci(2));
fprintf(fpCSV,cvsstr);

fprintf(fpCSV,'\n');
disp( sprintf('\n\n =======   All rats pooled:  S1B(CPP) and S2B(CPP) significance test:') );
tail = -1
[t, p, ci, ndof, h] = SlopeDiff_TTest(s2,sig2,n2,s4,sig4,n4, alpha, tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S1B(CPP) and S2B(CPP)',s2,s4,sig2,sig4,n2,n4,alpha,tail,t,p,h,s2-s4,ci(1),ci(2));
fprintf(fpCSV,cvsstr);
tail = -0
[t, p, ci, ndof, h] = SlopeDiff_TTest(s2,sig2,n2,s4,sig4,n4, alpha, tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S1B(CPP) and S2B(CPP)',s2,s4,sig2,sig4,n2,n4,alpha,tail,t,p,h,s2-s4,ci(1),ci(2));
fprintf(fpCSV,cvsstr);
tail = -0
[t, p, ci, ndof, h] = SlopeDiff_TTest(s2,sig2,n2,s4,sig4,n4, alpha, tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S1B(CPP) and S2B(CPP)',s2,s4,sig2,sig4,n2,n4,alpha,tail,t,p,h,s2-s4,ci(1),ci(2));
fprintf(fpCSV,cvsstr);

fprintf(fpCSV,'\n');
disp( sprintf('\n\n =======   All rats pooled:  S1B(saline) and S1B(CPP) significance test:') );
tail = -1
[t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s2,sig2,n2, alpha,tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S1B(saline) and S1B(CPP)',s1,s2,sig1,sig2,n1,n2,alpha,tail,t,p,h,s1-s2,ci(1),ci(2));
fprintf(fpCSV,cvsstr);
tail = +1
[t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s2,sig2,n2, alpha,tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S1B(saline) and S1B(CPP)',s1,s2,sig1,sig2,n1,n2,alpha,tail,t,p,h,s1-s2,ci(1),ci(2));
fprintf(fpCSV,cvsstr);
tail = 0
[t, p, ci, ndof, h] = SlopeDiff_TTest(s1,sig1,n1,s2,sig2,n2, alpha,tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S1B(saline) and S1B(CPP)',s1,s2,sig1,sig2,n1,n2,alpha,tail,t,p,h,s1-s2,ci(1),ci(2));
fprintf(fpCSV,cvsstr);


fprintf(fpCSV,'\n');
disp( sprintf('\n\n =======   All rats pooled:  S2B(saline) and S2B(CPP) significance test:') );
tail = -1
[t, p, ci, ndof, h] = SlopeDiff_TTest(s3,sig3,n3,s4,sig4,n4, alpha,tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S2B(saline) and S2B(CPP)',s3,s4,sig3,sig4,n3,n4,alpha,tail,t,p,h,s3-s4,ci(1),ci(2));
fprintf(fpCSV,cvsstr);
tail = +1
[t, p, ci, ndof, h] = SlopeDiff_TTest(s3,sig3,n3,s4,sig4,n4, alpha,tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S2B(saline) and S2B(CPP)',s3,s4,sig3,sig4,n3,n4,alpha,tail,t,p,h,s3-s4,ci(1),ci(2));
fprintf(fpCSV,cvsstr);
tail = 0
[t, p, ci, ndof, h] = SlopeDiff_TTest(s3,sig3,n3,s4,sig4,n4, alpha,tail)
cvsstr = sprintf('%s,%s,%9.4f,%9.4f,%9.4f,%9.4f,%d,%d,%9.4f,%d,%9.4f,%9.4f,%d,%9.4f,%9.4f,%9.4f,\n',...
    'All','S2B(saline) and S2B(CPP)',s3,s4,sig3,sig4,n3,n4,alpha,tail,t,p,h,s3-s4,ci(1),ci(2));
fprintf(fpCSV,cvsstr);




fclose(fpCSV);

diary off;
		
