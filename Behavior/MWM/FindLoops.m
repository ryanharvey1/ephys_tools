function [FinalLoop]=FindLoops(X,Y)
% FindLoops finds loops in XY data
%
% INPUT:
%           X & Y: column vectors containing x and y coordinants
%
% OUTPUT:
%           FinalLoop: Struct containing xy coordinants of loops
%
% by Ryan H., Laura B., Sep.2017

Z=[X,Y];
% FIND ALL POSSIBLE LOOPS
L=1;
for j=1:2
    for ii=round(1:10:length(Z)*2)
        if j==1; Z=[X,Y];else Z=[flipud(X),flipud(Y)];end
        for iframe=ii:length(Z)
            if iframe>length(Z); break; end
            for icom=iframe+1:length(Z)
                P=InterX([Z(iframe:icom,:)]');
                if ~isempty(P)==1 && length(Z(iframe:icom,1))>10                
                    start=iframe;
                    while isempty(P)==0 
                        P=InterX([Z(start:icom,:)]');
                        start=start+1;
                    end
                    if isempty(P)==1;start=start-2;end
                    LOOP=Z(start:icom,:);
                    if j==2; LOOP=[flipud(LOOP(:,1)),flipud(LOOP(:,2))];end
                    ALL_Loops(L).loops=LOOP; L=L+1;
                    Z(iframe:icom,:)=[];
                    break
                end
            end
        end
    end
end

if ~exist('ALL_Loops','var');FinalLoop=NaN;return;end

% FIND UNIQUE LENGTHS OF LOOPS
lengths=zeros(length(ALL_Loops),1);
for i=1:length(ALL_Loops)
    lengths(i,1)=length(ALL_Loops(i).loops);
end
[~,ia,~]=unique(lengths);

% STORE 
for i=1:length(ia)
    semiFinalLoop(i).loops=ALL_Loops(ia(i)).loops; 
end

% FIND UNIQUE LOOPS BY XY CORRELATION
for i=1:length(semiFinalLoop)
    for ii=i+1:length(semiFinalLoop)
        if isempty(semiFinalLoop(i).loops) || isempty(semiFinalLoop(ii).loops);continue;end
        sample=semiFinalLoop(i).loops;
        comp=semiFinalLoop(ii).loops;
        if size(sample,1)~=size(comp,1)
            if size(sample,1)-size(comp,1)>size(comp,1)-size(sample,1)
                comp=[comp;nan(size(sample,1)-size(comp,1),size(sample,2))];
            else
                sample=[sample;nan(size(comp,1)-size(sample,1),size(comp,2))];
            end
        end
        r=corr(sample,comp,'rows','complete');
        if sum(diag(r,0)>0.55)==2
            semiFinalLoop(ii).loops=[];
        end
    end
end

% REMOVE SHORT LOOPS
for i=1:length(semiFinalLoop)
    if isempty(semiFinalLoop(i).loops);continue;end
    if sum(sqrt(diff(semiFinalLoop(i).loops(:,1)).^2 + diff(semiFinalLoop(i).loops(:,2)).^2))<15
        semiFinalLoop(i).loops=[];
    end
end

% REMOVE EMPTY FIELDS
FinalLoop = semiFinalLoop(~cellfun(@isempty,{semiFinalLoop.loops}));






% \/ OLD CODE \/
% L=1;
% for j=1:2
%     for ii=round(1:4:length(Z)*2)
%         if j==1; Z=[X,Y];else; Z=[flipud(X),flipud(Y)];end
%         for iframe=ii:length(Z)
%             if iframe>length(Z); break; end
%             for icom=iframe+1:length(Z)
%                 if ismembertol(Z(iframe,:),Z(icom,:),0.085,'ByRows',true)==1 && length(Z(iframe:icom,1))>10
%                     LOOP=Z(iframe:icom,:);
%                     if j==2; LOOP=[flipud(LOOP(:,1)),flipud(LOOP(:,2))];end
%                     ALL_Loops(L).loops=LOOP; L=L+1;
%                     Z(iframe:icom,:)=[];
%                     break
%                 end
%             end
%         end
%     end
% end


% % FIND ALL LOOPS WITH 1 CROSSING
% for i=1:length(ALL_Loops)
%     P = InterX([ALL_Loops(i).loops]'); 
%     if isempty(P)==1
%         ALL_Loops(i).loops=[];
%         continue;
%     elseif size(P,2)~=1
%         ALL_Loops(i).loops=[];
%         continue;
%     end 
% end
% ALL_Loops=ALL_Loops(~cellfun(@isempty,{ALL_Loops.loops}));
% if isempty(ALL_Loops);FinalLoop=NaN;return;end
