function [coherenceout]=tiltedcoherence(ratemap)
% Coherence

[r,c]=size(ratemap);

ratemap=[NaN(1,c);ratemap;NaN(1,c)];

[r,c]=size(ratemap);


ratemap=[NaN(r,1),ratemap,NaN(r,1)];

rawpixel=[];
avgneighbor=[];

for c=2:size(ratemap,2)-1
    for r=2:size(ratemap,1)-1
        % extract raw pixel
        rawpixel=[rawpixel;ratemap(r,c)];
        % find 8 neighbors
        neighbors(1) = ratemap(r-1,c-1); % Upper left.  r = row, c = column.
        neighbors(2) = ratemap(r-1,c); % Upper middle.  r = row, c = column.
        neighbors(3) = ratemap(r-1,c+1); % Upper right.  r = row, c = column.
        neighbors(4) = ratemap(r,c-1); % left.  r = row, c = column.
        neighbors(5) = ratemap(r,c+1); % right. r = row, c = column.
        neighbors(6) = ratemap(r+1,c+1); % Lowerleft.  r = row, c = column.
        neighbors(7) = ratemap(r+1,c); % lower middle.  r = row, c = column.
        neighbors(8) = ratemap(r+1,c-1); % Lower left.  r = row, c = column.
        % average them
        avgneighbor=[avgneighbor;nanmean(neighbors)];
    end
end

coherenceout=corr(rawpixel,avgneighbor);
end


