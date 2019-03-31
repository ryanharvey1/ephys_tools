function [C,B] = SpectCrossCorr(S,t,nbBins)

if round(nbBins/2) == nbBins/2
	nbBins = nbBins+1;
end

dS = Data(S)';
tS = Range(S);
totBin = size(dS,2);

C = zeros(size(dS,1),nbBins);
dt = median(diff(Range(S)));
l = (nbBins-1)/2;


%  try
for i=1:length(t)

	trigT = binsearch_floor(tS,t(i));
	if trigT>l & trigT<totBin-l
		C = C+dS(:,trigT-l:trigT+l);
	else
%  		warning('size error')
%  		keyboard
	end

end
%  catch
%  keyboard
%  end
C = C/length(t);
B = [-l*dt:dt:l*dt];