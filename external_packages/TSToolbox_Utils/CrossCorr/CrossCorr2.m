function [C,B] = CrossCorr(t1,t2,binsize,nbins)

nt1  = length(t1);
nt2 = length(t2);

% we want nbins to be odd */
if floor(nbins / 2)*2 == nbins
    nbins = nbins+1;
end
  
m = - binsize * ((nbins+1) / 2);
B = zeros(nbins,1);
for j = 1:nbins
	B(j) = m + j * binsize;
end
    
% cross correlations */
  
w = ((nbins) / 2) * binsize;
C = zeros(nbins,1);
i2 = 2;

for i1 = 1:nt1
    lbound = t1(i1) - w;
    while t2(i2) < lbound && i2 < nt2
        i2 = i2+1;
    end
    while t2(i2-1) > lbound && i2 > 2
	i2 = i2-1;
    end
    
    rbound = lbound;
    l = i2;
    for j = 1:nbins
        k = 0;
        rbound = rbound+binsize;
        while t2(l) < rbound && l < nt2-1  
            l = l+1;
            k = k+1;
        end
      
	  C(j) = C(j)+k;
    end
  
end


  for j = 1:nbins
    C(j) = C(j)/(nt1 * binsize / 10000);
 
  end
  
