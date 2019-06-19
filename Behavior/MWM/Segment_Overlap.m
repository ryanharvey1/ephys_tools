function segments=Segment_Overlap(X,Y,lseg,ovlp)
%Segment_Overlap Splits the trajectory in segments of length
% lseg with an overlap of ovlp %
% Returns an structure of xy coords for each segment
% 
% From MWM-ML-GEN, edited by Ryan H 2017

% figure(1);plot(X,Y);hold on

Z=[X,Y];
n = size(Z,1);

% compute cumulative distance vector
cumdist = zeros(1, n);
for i = 2:n
    cumdist(i) = nansum([cumdist(i - 1); norm(Z(i,:) - Z(i-1,:) )]); %changed to nansum lb 06/2019 
    % nansum allows the vector to contain NaN. May be problematic for accurate segmentation with a large amount of contiguous nans.                                                                  
    %                                   
end

% step size
off = lseg*(1. - ovlp);
% total number of segments - at least 1
if cumdist(end) > lseg
    nseg = ceil((cumdist(end) - lseg) / off) + 1;
    off = off + (cumdist(end) - lseg - off*(nseg - 1))/nseg;
else
    nseg = 1;
end
% segments are trajectories again -> construct empty object
%     segments = trajectories([]);

for seg = 0:(nseg - 1)
    starti = 0;
    seg_off = 0;
    pts = [];
    if nseg == 1
        % special case: only 1 segment, don't discard it
        pts = Z;
    else
        for i = 1:n
            if cumdist(i) >= seg*off
                if starti == 0
                    starti = i;
                end
                if cumdist(i) > (seg*off + lseg)
                    % done we are
                    break;
                end
                if isempty(pts)
                    seg_off = cumdist(i);
                end
                % otherwise append point to segment
                pts = [pts; Z(i, :)];
            end
        end
    end
    segments(seg+1).segments=pts;
%     figure(1);plot(pts(:,1),pts(:,2),'linewidth',3);pause(1)
end
end