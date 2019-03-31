function [U2,p] = watson1962(N1, N2)
%G.S. Watson's statistic for the Goodness-of-fit tests on a circle, pt
%II.
%See
%http://biomet.oxfordjournals.org/content/49/1-2/57.full.pdf
%for reference.

%Sample Data:
% N1 = [50; 290; 300; 300; 305; 320; 330; 330; 335; 340; 340; 355];
% N2 = [70; 155; 190; 195; 215; 235; 235; 240; 255; 260; 290; 300; 300; 300];


N1 = sort(N1);
N2 = sort(N2);

N1_set = zeros(size(unique(N1),1),1);
N2_set = zeros(size(unique(N2),1),1);

N1_occurs = zeros(size(unique(N1),1),1);
N2_occurs = zeros(size(unique(N2),1),1);

for i = 1:size(N1)
    N1_occurs(i) = sum(N1==N1(i));    
end

N1_set = unique([N1 N1_occurs],'rows');

for i = 1:size(N2)
    N2_occurs(i) = sum(N2==N2(i));    
end

N2_set = unique([N2 N2_occurs],'rows');

all_N = sort(unique([N1_set(:,1); N2_set(:,1)]));

N1 = size(N1,1);
N2 = size(N2,1);

N1_count = 0;
N2_count = 0;

z_values = zeros(size(all_N,1),1);
N_k = zeros(size(all_N,1),1);
% Note: N1_set is unique values. 

for i = 1:size(all_N,1)
    el = all_N(i);
    N1_idx = find(N1_set(:,1) == el);
    N2_idx = find(N2_set(:,1) == el);
    
    if  ~isempty(N1_idx) && ~isempty(N2_idx)
        N1_count = N1_count + N1_set(N1_idx,2);
        N2_count = N2_count + N2_set(N2_idx,2);
        N_k(i) = N1_set(N1_idx,2) + N2_set(N2_idx,2);
    elseif ~isempty(N1_idx)
        N1_count = N1_count + N1_set(N1_idx,2);
        N_k(i) = N1_set(N1_idx,2);
    elseif ~isempty(N2_idx)
        N2_count = N2_count + N2_set(N2_idx,2);
        N_k(i) = N2_set(N2_idx,2);
    else
        disp('something went wrong');
    end
    z = (N1_count / N1) - (N2_count / N2);
    z_values(i) = z;
end

S1 = sum(N_k .* z_values);
S2 = sum(N_k .* z_values.^2);

U2 = ((N1*N2) / ((N1+N2)^2)) * (S2 - (S1^2)/(N1+N2));
p = 2*exp(19.74*-U2);
end