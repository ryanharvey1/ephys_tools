function H = hsvOccMap(occH,varargin)

%Adrien Peyrache, 2012

slope = 10;

offSet = 0.7;
if ~isempty(varargin)
    offSet = varargin{1};
end
H = zeros(size(occH,1),size(occH,2),3);

%pure empirical definition that works best
if min(occH) == max(occH)
    occH = ones(size(occH));
else
    
    occH = (occH-min(occH(:)))/(max(occH(:))-min(occH(:)));
    occH = tanh(slope*occH)+offSet;
    occH = occH/max(occH(:));
end

H(:,:,3) = occH;
H(:,:,2) = 1;
