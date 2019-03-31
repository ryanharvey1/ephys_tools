function [HalfWidth,field]=FindHalfWidth(HDtuning)
%by Ryan H. edits by Laura B. 16-Mar-18

if sum(HDtuning)==0
    HalfWidth=0;
    field=zeros(1,length(HDtuning));
    return
end

cm2removeR=0;cm2removeL=0;

[M,I]=max(HDtuning); Thres=M*0.50;
middlebin=round(median(1:length(HDtuning)));
HDtuning=circshift(HDtuning,(middlebin-I)-1);

[M,I]=max(HDtuning);
field=zeros(1,length(HDtuning)); field(I)=1;



% FORWARD
for i=1:length(HDtuning)-I
    if HDtuning(I)~=length(HDtuning)
        if HDtuning(I+i)>Thres
            field(I+i)=1;
        elseif HDtuning(I+i)<Thres
            cm2removeR=(length(HDtuning)-(I+(i-1)))*5;
            break
        end
    end
end

% BACKWARD
for i=1:I-1
    if HDtuning(I)~=1
        if HDtuning(I-i)>Thres
            field(I-i)=1;
        elseif HDtuning(I-i)<Thres
            cm2removeL=((I-(i-1))-1)*5;
            break
        end
    end
end

% frames=[rescale(matrix(:,2),1,length(HDtuning)),rescale(matrix(:,3),1,2)];
% 
% [row,col]=find([field;field]);k=boundary(row,col);bound=[col(k),row(k)];
% in=inpolygon(frames(:,1),frames(:,2),bound(:,1),bound(:,2));
% spks_VEL=matrix(logical(in) & matrix(:,6)==1,:);
% 
% 
% %
% %                 pix2removeR=cm2removeR/pixelDist;
% %                 pix2removeL=cm2removeL/pixelDist;
% %
% %                 % FIND MAX AND MIN IN X DIM
% %                 x_min=min(matrix(:,2)); x_max=max(matrix(:,2));
% %                 % FIND VALUE OF PIXEL CUT OFF
% %                 lowercut=x_min+pix2removeL; highercut=x_max-pix2removeR;
% %                 % CUT OUT PIXELS
% %                 matrix=matrix((matrix(:,2)>lowercut & matrix(:,2)<highercut),:);
% %
% %                 spks_VEL = matrix(matrix(:,6) == 1,:);
HalfWidth=sum(field)*(360/length(HDtuning));
end