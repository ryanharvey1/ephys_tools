function [x,y] = InterpolateTrackerZeros_v2(x,y)

i0 = find(x+y == 0);
if isempty(i0)
    % no corrections necessary
    return;
end%if

ndata = length(x);

% check if last point is zero; if yes, then truncate all trailing zeros!
while(x(ndata) == 0 & y(ndata) == 0)
    x(ndata)=[];
    y(ndata)=[];
    ndata = length(x);
end%if

% interpolate linearly each block of consecutive zeros
for i = 1:ndata
    if(x(i) == 0 & y(i) == 0)
        nzeros = 1;
        ibstart = i-1;
        i=i+1;
        while (x(i) == 0 & y(i) == 0)
            % count consecutive zeros
            nzeros = nzeros + 1;
            i=i+1;
        end%while
        ibend = i;
        x(ibstart:ibend)=linspace(x(ibstart),x(ibend),nzeros+2);
        y(ibstart:ibend)=linspace(y(ibstart),y(ibend),nzeros+2);
    end%if
end%for
end