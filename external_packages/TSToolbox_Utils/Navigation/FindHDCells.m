function FindHDCells(varargin)

% This function computes some Head Direction cell statistics
% Must be launched from the folder containing data
% Folder must have the same name as position file basename.

% Adrien Peyrache 2012

epIx = [];
if ~isempty(varargin)
    epIx = varargin{1};
end

if length(varargin)==2
    fname = varargin{2};
    fname = ['_' fname];
else
    fname = '';
end

[~,fbasename,~] = fileparts(pwd);
B = [];
AngHisto = [];
% 
% if exist('Analysis/BehavEpochs.mat','file')
%     load('Analysis/BehavEpochs.mat');
% end

load('Analysis/SpikeData.mat','S','shank');
nbC = length(S);

if ~isempty(epIx)   
    if isa(epIx,'char')
        load('Analysis/BehavEpochs.mat',epIx);
        eval(['wakeEp = ' epIx]);
        if isempty(fname)
            fname = ['_' epIx];
        end
    else
        if epIx >0
            wakeEp = LoadEpochFromDatSeg(epIx);
        else
            wakeEp = timeSpan(S);
        end
    end
else
    load('Analysis/BehavEpochs.mat');
end

hdCellStats = NaN(nbC,5);

if exist('wakeEp','var')

[ang,angGR] = HeadDirection_Wrapper(fbasename);
if isempty(epIx)
    [X,Y,~,wstruct] = LoadPosition_Wrapper(fbasename);
else
    [X,Y,~,wstruct] = LoadPosition_Env(fbasename,epIx);
end

[linSpd,angSpd] = LoadSpeed_Wrapper(fbasename,wstruct); 
%calling wstruct prevents reading again the whl file

wakeEp = intersect(wakeEp,angGR);

mobilityV = thresholdIntervals(linSpd,2,'Direction','Above');
wakeEp = intersect(wakeEp,mobilityV);
epAng = intersect(angGR,wakeEp);

ang = Restrict(ang,epAng);
S = Restrict(S,wakeEp);

meanAng     = zeros(nbC,1);
multiPks    = zeros(nbC,1);
kappa       = zeros(nbC,1);
pVal        = zeros(nbC,1);
AngHisto    = [];
peakFr      = zeros(length(S),1);
hdStability = zeros(length(S),1);

for ii=1:nbC
    if any(Range(Restrict(S{ii},epAng)))
        [AngHisto(:,end+1),B,meanAng(ii),pVal(ii),kappa(ii)] = HeadDirectionField_Norm(S{ii},ang,epAng);
        
        pks = LocalMinima(-AngHisto(:,end),60,-1);
        [peakFr(ii),mxIx] = max(AngHisto(:,end)); 
%         watsonU(ii) = watsons_U2(Data(Restrict(ang,S{ii})),Data(ang));
%         watsonU(ii)

        epHalf          = regIntervals(epAng,2);
        h1              = HeadDirectionField_Norm(S{ii},ang,epHalf{1});
        h2              = HeadDirectionField_Norm(S{ii},ang,epHalf{2});        
        hdS             = corrcoef(h1,h2);
        hdStability(ii)  = hdS(1,2);

        if length(pks)>1
            if sum(AngHisto(pks,end)/peakFr(ii) > 0.25)>1
                multiPks(ii) = 1;
            end
        end

        if 0
            figure(1),clf
            polar(B,AngHisto(:,end))
            disp([meanAng(ii) pVal(ii) kappa(ii)])
            pause
        end
    end
end

meanAng = mod(meanAng,2*pi);
hdIndex = [];
phdIndex = [];
hdScore = [];

[hdIndex,phdIndex] = HDIndex(S,ang,epAng);
%hdScore = HDScore(AngHisto,meanAng,B(1:end-1));

%hdIx = pVal<0.001 & kappa>1 & peakFr>1;

hdIx = peakFr>1 & phdIndex<=0.001 & kappa>1 & hdStability>0.75;
disp(['Total # of HD cells:' num2str(sum(hdIx))])
hdCellStats = [meanAng phdIndex kappa peakFr hdIx];


else
    keyboard
    %hdCellStats = NaN(nbC,5);
end

moveEp = wakeEp;
info = {'No Info';'No Info';'1st column, preferred head direction; 2nd, p-value; 3rd, kappa value; 4th: peak f.r. in field';'moving epoch';'multi peaks?';'hdIndex';'phdIndex';'hdStability'};
SaveAnalysis(pwd,['HDCells' fname],{AngHisto;B;hdCellStats;moveEp;multiPks;hdIndex;phdIndex;hdStability},{'AngHisto';'B';'hdCellStats';'moveEp';'multiPks';'hdIndex';'phdIndex';'hdStability'},info);