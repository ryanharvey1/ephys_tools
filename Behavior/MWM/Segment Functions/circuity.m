function [c]=circuity(x,y)
obs=sum(sqrt(diff(x).^2+diff(y).^2));
TL=sqrt((x(1,1)-x(end,1))^2+(y(1,1)-y(end,1)).^2);
c=TL/obs;
end