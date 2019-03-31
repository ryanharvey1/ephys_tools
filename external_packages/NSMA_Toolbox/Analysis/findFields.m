function PF = findFields(Qc, bins, plotflag, minFR, troughtopeak, propFRlimit)
% Find Place Fields (using only smoothed Firing Rate tuning curve)
% PF = findFields(Qc, binsize, plotflag, minFR, troughtopeak, propFRlimit)
%
% inputs:    Qc: SMOOTHED matrix of place field tuning curves (cells in
%                rows, position bins in columns)
%            bins: vector of the same lenght as a row of Qc, detefining the
%                  position value of each bin (used for plotting only)
%            plotflag (optional): if 0, don't plot anything, if 1 (default) plot each
%                                 cell, if 2 plot each place field (not
%                                 recommended for many cells!)
%            other optional variables to effect how fields are defined, see below
%
% Place field defining algorithm:
%   manipulatible (optional) variables:
%       minFR   (0.8Hz)     estimate of baseline firing rate
%       troughtopeak (0.7)  ratio of trough to peak at which to divide a
%                           field into two
%       propFRlimit (0.1)   proportion of max firing rate of field that can
%                           get cut off
%   use firing rate binned to position along track and SMOOTHED
%   find peaks > minFR and separated by at least minFR
%   go through each peak to define its boundaries:
%       find troughs out from 25% of peak
%       pick the first trough away from peak:
%           following which the firing rate is less than minFR OR propFRlimit*peakFR
%           OR peaktotrough*nextpeak is greater than the trough
%   for any place fields that share a boundary, check if trough:peak <
%       trough to peak, and if so merge them
%
% output:
%   PF - cell array including a matrix for each cell, which includes a row
%   for each place field, with the first column containing the start
%   position of each field, the second column the end position of each
%   field, and the third column the peak FR of each field. 
%
% uses peakdetz.m to find peaks and troughs
%
% ZN 04/2010

[nCells, nBins]=size(Qc);

if nargin<6
    propFRlimit=0.1;      % proportion of max firing rate of field that can get cut off
    if nargin<5
        troughtopeak=0.5;      % ratio of trough to peak at which to divide a field into two
        if nargin<4
            minFR=0.8;      % estimate of baseline firing rate. used as cut-off for finding firing rate peaks and troughs for boundaries
            if nargin<3
                plotflag=1;
                if nargin<2
                    bins=1:nBins;
                end
            end
        end
    end
end

if plotflag && length(bins)~=nBins
    disp('Warning: bins entered do not match bins in Qc')
    bins=1:nBins;
end

if plotflag
    plotcell=1;
    if plotflag==2
        plotPF=1;
    else
        plotPF=0;
    end
else
    plotPF=0;
    plotcell=0;
end

PF=cell(nCells,1);      % cell array for storing place field info

for c=1:nCells
    if plotcell
        cf=figure;
        plot(bins, Qc(c,:), 'k-')
        axis tight
        ylabel('Normalized firing rate')
    end
    
    % find firing rate peaks
    [peak,trough]=peakdetz(Qc(c,:), minFR, 1);
    
    if plotcell
        if ~isempty(peak)
            hold on; plot(bins(peak(:,1)), Qc(c,peak(:,1)),'r.')
        end
        if ~isempty(trough)
            hold on; plot(bins(trough(:,1)), Qc(c,trough(:,1)),'b.')
        end
    end
    
    [npeaks,w]=size(peak);
    [ntroughs,w]=size(trough);
    
    % find place field boundaries for each peak
    PF{c}=zeros(npeaks,3);     % matrix for storing start (1) and end (2) position and peaks (3) of each PF
    for p=1:npeaks
        ipeak=peak(p,1);
        peakFR=peak(p,2);
        if p<2
            ileftbound=1;
        else
            ileftbound=trough(p-1,1);
        end
        
        if p>ntroughs
            irightbound=nBins;
        else
            irightbound=trough(p,1);
        end
        
        if plotPF
            pf=figure; plot(ileftbound:irightbound, Qc(c,ileftbound:irightbound), 'k')
        end
        
        % find boundaries of 25% of peak FR of field
        ileftside=ileftbound-1 + find(Qc(c,ileftbound:ipeak)<0.25*peakFR, 1, 'last');
        irightside=ipeak-1 + find(Qc(c,ipeak:irightbound)<0.25*peakFR, 1, 'first');
        
        if plotPF
            hold on; plot([ileftside irightside], [Qc(c,ileftside) Qc(c,irightside)], 'b.')
        end
        
        % find small troughs on left side of field
        if isempty(ileftside)
            [p2,t2]=peakdetz(Qc(c,ileftbound:ipeak), 0.01, 0, 1);
            %disp(num2str(p));
            if ~isempty(t2)
                t2(:,1)=ileftbound-1 + t2(:,1);
                if plotPF, hold on; plot(t2(:,1), t2(:,2), 'g.'); end
            end
        else
        	[p2,t2]=peakdetz(Qc(c,ileftbound:ileftside), 0.01, 0, 1);
            if ~isempty(t2)
                t2(:,1)=ileftbound-1 + t2(:,1);
                if plotPF, hold on; plot(t2(:,1), t2(:,2), 'g.'); end
            end
        end
        
        % go through small troughs and find the one after which firing
        % rate<minFR (or propFRlimit*peakFR), or troughtopeak* next max is 
        % greater than the trough
        [nummin,w]=size(t2);
        [nummax,w]=size(p2);
        newmin=1;
        if nummin>0;
            inewmin=t2(newmin,1);
        end
        while newmin<=nummin && (max(Qc(c,ileftbound:inewmin))>minFR || max(Qc(c,ileftbound:inewmin))>propFRlimit*peakFR)
            if newmin>nummax
                nextpeak=0;
            else
                nextpeak=p2(newmin,2);
            end
            if nextpeak*troughtopeak>t2(newmin,2) || (nextpeak<propFRlimit*peakFR && nextpeak<minFR)
                break;
            end
            if plotPF, hold on; plot(t2(newmin,1), t2(newmin,2), 'r.'); end
            newmin=newmin+1;
            if newmin<=nummin
                inewmin=t2(newmin,1);
            end
        end
        % set new left limit
        if newmin<=nummin
            ileftbound=t2(newmin,1);
        end
        
        if plotPF, hold on; plot([ileftbound ileftbound], [0 Qc(c,ipeak)], 'r-'); end
        
        % find small troughs on right side of field
        if isempty(irightside)
            [p2,t2]=peakdetz(Qc(c,ipeak:irightbound), 0.01, 0, 0);
            %disp(num2str(p));
            if ~isempty(t2)
                t2(:,1)=ipeak-1 + t2(:,1);
                if plotPF, hold on; plot(t2(:,1), t2(:,2), 'g.'); end
            end
        else
        	[p2,t2]=peakdetz(Qc(c,irightside:irightbound), 0.01, 0, 0);
            if~isempty(t2)
                t2(:,1)=irightside-1 + t2(:,1);
                if plotPF, hold on; plot(t2(:,1), t2(:,2), 'g.'); end
            end
        end
        
        % go through small troughs and find the one after which firing
        % rate<minFR (or propFRlimit*peakFR), or troughtopeak* next max is 
        % greater than the trough
        [nummin,w]=size(t2);
        [nummax,w]=size(p2);
        newmin=1;
        if nummin>0;
            inewmin=t2(newmin,1);
        end
        while newmin<=nummin && (max(Qc(c,inewmin:irightbound))>minFR || max(Qc(c,ileftbound:inewmin))>propFRlimit*peakFR)
            if plotPF, hold on; plot(t2(newmin,1), t2(newmin,2), 'r.'); end
            if newmin>nummax
                nextpeak=0;
            else
                nextpeak=p2(newmin,2);
            end
            if nextpeak*troughtopeak>t2(newmin,2) || (nextpeak<propFRlimit*peakFR && nextpeak<minFR)
                break;
            end
            newmin=newmin+1;
            if newmin<=nummin
                inewmin=t2(newmin,1);
            end
        end
        % set new right limit
        if newmin<=nummin
            irightbound=t2(newmin,1);
        end
        
        if plotPF
            hold on; plot([irightbound irightbound], [0 Qc(c,ipeak)], 'r-')
            hold off;
        end
        
        PF{c}(p,:)=[ileftbound, irightbound, peakFR];
    end
    
    if plotcell && ~isempty(PF{c})
        figure(cf);
        hold on; plot([bins(PF{c}(:,1)); bins(PF{c}(:,1))], [0 max(Qc(c,:))], 'g--')
        plot([bins(PF{c}(:,2)); bins(PF{c}(:,2))], [0 max(Qc(c,:))], 'g-')
    end
    
    % check fields that are touching and merge them if appear as 1 field
    pf=1;
    while pf<npeaks
        if PF{c}(pf,2)==PF{c}(pf+1,1)
            % check if trough between fields is less than troughtopeak* 
            % smallest peak. if so, merge the two fields
            if Qc(c,PF{c}(pf,2))>troughtopeak*min(PF{c}(pf:pf+1,3))
                PF{c}(pf,2)=PF{c}(pf+1,2);
                PF{c}(pf,3)=max(PF{c}(pf:pf+1,3));
                PF{c}(pf+1,:)=[];
                npeaks=npeaks-1;
            end
        end
        pf=pf+1;
    end
    
    if plotcell && ~isempty(PF{c})
        figure(cf);
        hold on; plot([bins(PF{c}(:,1)); bins(PF{c}(:,1))], [0 max(Qc(c,:))], 'r--')
        plot([bins(PF{c}(:,2)); bins(PF{c}(:,2))], [0 max(Qc(c,:))], 'r-')
    end
end

