function [D,M]=Displacement(Ratemap1,Ratemap)
% Displacement returns the field displacement between two ratemaps in
% degrees of rotation & rotation correlation
% 
% INPUT:    Ratemap1    Sample ratemap. For example, session 1 ratemap before cue shift
%           Ratemap     Comparison ratemap. For example, session 2 ratemap after cue shift
% 
% OUTPUT:   D           Degree of rotation in field
%           M           Correlation to sample ratemap given rotation
% 
% RYAN HARVEY 7/24/17

rALL=[];
% PAD MAPS
Ratemap1=padarray(Ratemap1,[2 2]);
Ratemap=padarray(Ratemap,[2 2]);

% FORCE MAPS TO BE SAME SIZE
if size(Ratemap1,1)~=size(Ratemap,1)
    if size(Ratemap1,1)-size(Ratemap,1)>size(Ratemap,1)-size(Ratemap1,1)
        Ratemap=[Ratemap;zeros(size(Ratemap1,1)-size(Ratemap,1),size(Ratemap,2))]; 
    else
        Ratemap1=[Ratemap1;zeros(size(Ratemap,1)-size(Ratemap1,1),size(Ratemap1,2))];
    end
end
if size(Ratemap1,2)~=size(Ratemap,2)
    if size(Ratemap1,2)-size(Ratemap,2)>size(Ratemap,2)-size(Ratemap1,2)
        Ratemap=[Ratemap,[zeros(size(Ratemap1,2)-size(Ratemap,2),size(Ratemap,1))]'];
    else
        Ratemap1=[Ratemap1,[zeros(size(Ratemap,2)-size(Ratemap1,2),size(Ratemap1,1))]'];
    end
end
% ROTATE MAP IN 6 DEGREE BINS & CORRELATE TO SAMPLE MAP
[rows,cols]=size(Ratemap);
c=[round(size(Ratemap)/2)]';
[i,j]=find(Ratemap);I=[i j]';
c=c*ones(1,length(i));
for angle=0:6:359
    theta = deg2rad(angle);
    R=[cos(theta) -sin(theta);sin(theta) cos(theta)];
    I2=round(R*(I-c)+c);
    if sum(I2(1,:)>rows) || sum(I2(2,:)>cols)
        for ii=1:length(I2(1,:))
            if I2(1,ii)>rows;I2(1,ii)=I2(1,ii)-(I2(1,ii)-rows);end
            if I2(2,ii)>cols;I2(2,ii)=I2(2,ii)-(I2(2,ii)-cols);end
        end
    end
    I2(I2==0)=1; I2=abs(I2);
    r=corr2(Ratemap1,full(sparse(I2(1,:),I2(2,:),1,rows,cols))); rALL = [rALL; r];
end
% FIND MAX CORRELATION & DEGREE SHIFT
[M,Index]=max(rALL);D=Index*6;
end