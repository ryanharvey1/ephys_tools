%  Pennartz_RegAnalysis script
% 
%  Loads in the values from ABM (all bias matrix)
%  and runs the matlab function regress on the S1-B
%  and S2-B bias values.
%  This returns the r2 and p values which are saved.
group = 'CPP'    % set to 'Saline' for saline, 'CPP' for CPP
%group = 'Saline'    % set to 'Saline' for saline, 'CPP' for CPP
runNrStr = 'Run000_'
outDir = 'D:\CPP1\JasonCPPSleepStudy\jason\Out'

plotRats = 0;
showErrorBars = 1;

if strcmpi(group,'Saline')
    run_label = 'Saline';
    indir = [outDir '\TempBias_Saline_15min'];
    ratlist = {'6318_200ms';'6401_200ms';'6495_200ms';'6496_200ms';'6568_200ms';'6591_200ms'};
elseif strcmpi(group,'CPP')   
    run_label = 'CPP';
    indir = [outDir '\TempBias_CPP_15min'];
    ratlist = {'6318_200ms';'6401_200ms';'6481_200ms';'6495_200ms';'6496_200ms';'6568_200ms';'6591_200ms'};
end


outputDir = indir
outPrefix = [runNrStr 'RegAnalysis_'];

diary off;                  % in case it is on from a bombed previous run
diary([outputDir, filesep, outPrefix, 'run_', strrep(datestr(datenum(now),0),':','-'),'.txt']);

%for diary output
run_label
group
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
        if showErrorBars
            hold on
            [P,stat] = polyfit(X(:,2),ABM(:,1),1);
            xx = [-1:0.2:1];
            [yy,delta] = polyval(P,xx,stat);
            errorbar(xx,yy,delta,'r');
            hold off
        end
        grid on
        [NDots,Ncols] = size(ABM)
        title([runNrStr ' ' ratlist{jj} ': B vs S1    (' num2str(NDots) ' data points)'])
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
        if showErrorBars
            hold on
            [P,stat] = polyfit(X(:,2),ABM(:,3),1);
            xx = [-1:0.2:1];
            [yy,delta] = polyval(P,xx,stat,'r');
            errorbar(xx,yy,delta);
            hold off
        end
        grid on
        title([runNrStr ' ' ratlist{jj} ': B vs S2    (' num2str(NDots) ' data points)'])
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
    
end

fprintf('======================================================================\n');
fprintf('======================================================================\n');
fprintf('============ All Rats All Session pooled together :===================\n\n');
[NDots,Ncols] = size(S1_bias_pooled);
NDots 

X = [ones(size(B_bias_pooled,1),1) B_bias_pooled];
[S1Bstats.b,S1Bstats.bint,S1Bstats.r,S1Bstats.rint,S1Bstats.stats] = regress(S1_bias_pooled,X);  % S1-B 
disp('S1Bstats.b ='); S1Bstats.b
disp('S1Bstats.bint ='); S1Bstats.bint
disp('S1Bstats.stats ='); S1Bstats.stats

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
if showErrorBars
    hold on
    [P,stat] = polyfit(X(:,2),S1_bias_pooled,1);
    xx = [-1:0.2:1];
    [yy,delta] = polyval(P,xx,stat);
    errorbar(xx,yy,delta,'r');
    hold off
end
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
if showErrorBars
    hold on
    [P,stat] = polyfit(X(:,2),S2_bias_pooled,1);
    xx = [-1:0.2:1];
    [yy,delta] = polyval(P,xx,stat);
    errorbar(xx,yy,delta,'r');
    hold off
end
grid on
title([runNrStr ' ' run_label ': ' num2str(length(ratlist)) ' Rats combined: B vs S2    (' num2str(NDots) ' data points)'])
xlabel('Behavior')
ylabel('Sleep 2')
axis([-1,1,-1,1])
axis equal

saveas(gcf,[outputDir, filesep, outPrefix, run_label,'_AllRatsCombined'],'fig');

fprintf('================= Done ================================================\n');

diary off;

