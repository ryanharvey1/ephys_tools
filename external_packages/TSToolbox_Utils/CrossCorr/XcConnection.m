function [SynC,ConStr] = XcConnection(HO,Hm,HeI,HeS,Hstd,B,varargin)

%  [H0,Hm,HeI,HeS,P] = XcConnection(t1,t2,binS,nbBins)
% Find Synpatic connectivity based on cross-correlograms
%  
%  INPUT:
%    H0: raw CrossCorr
%    HeI: lower bound of confidence Interval
%    HeS: upper bound  bound of confidence Interval
%  
%  OUTPUT:
%    Exc: 0 if no ecitation, 1 if cell1 excites cell2, 2 if cell2 excites cell1, 3 if both
%    Inh: 0 if no inhibition, 1 if cell1 inhibites cell2, 2 if cell2 inhibites cell1, 3 if both
%  
%  
%  Adrien Peyrache, based on Fujisawa et al. 2008

SynC = [0 0];
ConStr = [0 0];

if length(varargin)
  intC = varargin{1};
  minB = intC(1);
  maxB = intC(2);
else
  minB = 1;
  maxB = 4;
end

bIx = find(B>=minB & B<=maxB);
b = B(bIx);

h = HO(bIx);
m = Hm(bIx);
s = Hstd(bIx);
[v,ix] = max(h-HeS(bIx));
if v>0
  SynC(1) = 1;
  ConStr(1) = (h(ix)-m(ix))/s(ix);
end

[v,ix] = min(h-HeI(bIx));
if v<0
  p = (m(ix)-h(ix))/s(ix);
  if ~SynC(1) | p<ConStr(1)
      SynC(1) = -1;
      ConStr(1) = p;
  end
end


bIx = find(B>=-maxB & B<=-minB);
b = B(bIx);

h = HO(bIx);
m = Hm(bIx);
s = Hstd(bIx);
[v,ix] = max(h-HeS(bIx));
if v>0
  SynC(2) = 1;
  ConStr(2) = (h(ix)-m(ix))/s(ix);
end

[v,ix] = min(h-HeI(bIx));
if v<0
  p = (m(ix)-h(ix))/s(ix);
  if ~SynC(2) | p<ConStr(2)
      SynC(2) = -1;
      ConStr(2) = p;
  end
end
