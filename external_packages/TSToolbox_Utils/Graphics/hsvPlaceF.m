function pfh = hsvPlaceF(pf,varargin)

if ~isempty(varargin)
    m = varargin{1};
    mn = m(1);
    mx = m(2);
else
    mn = min(pf(:));
    mx = max(pf(:));
end
% pf = gausssmooth(pf,20);
pfh = 1-((pf-mn)/(1.5*(mx-mn))+0.333);