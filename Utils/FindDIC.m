function [spks_VEL,FieldWidth,field]=FindFF(Ratemap,matrix,pixelDist,track_length)

if sum(Ratemap)==0
    spks_VEL=zeros(1,size(matrix,2));
    FieldWidth=0;
    field=zeros(1,length(Ratemap));
    return
end

cm2removeR=0;cm2removeL=0;

[M,I]=max(Ratemap); Thres=M*0.20;
field=zeros(1,length(Ratemap)); field(I)=1;

% FORWARD
for i=1:length(Ratemap)-I
    if Ratemap(I)~=length(Ratemap)
        if Ratemap(I+i)>Thres;
            field(I+i)=1;
        elseif Ratemap(I+i)<Thres;
            cm2removeR=(length(Ratemap)-(I+(i-1)))*5;
            break
        end
    end
end

% BACKWARD
for i=1:I-1
    if Ratemap(I)~=1
        if Ratemap(I-i)>Thres;
            field(I-i)=1;
        elseif Ratemap(I-i)<Thres;
            cm2removeL=((I-(i-1))-1)*5;
            break
        end
    end
end

frames=[rescale(matrix(:,2),1,length(Ratemap)),rescale(matrix(:,3),1,2)];

[row,col]=find([field;field]);k=boundary(row,col);bound=[col(k),row(k)];
in=inpolygon(frames(:,1),frames(:,2),bound(:,1),bound(:,2));
spks_VEL=matrix(logical(in) & matrix(:,6)==1,:);


%
%                 pix2removeR=cm2removeR/pixelDist;
%                 pix2removeL=cm2removeL/pixelDist;
%
%                 % FIND MAX AND MIN IN X DIM
%                 x_min=min(matrix(:,2)); x_max=max(matrix(:,2));
%                 % FIND VALUE OF PIXEL CUT OFF
%                 lowercut=x_min+pix2removeL; highercut=x_max-pix2removeR;
%                 % CUT OUT PIXELS
%                 matrix=matrix((matrix(:,2)>lowercut & matrix(:,2)<highercut),:);
%
%                 spks_VEL = matrix(matrix(:,6) == 1,:);
FieldWidth=sum(field)*(track_length/length(field));
end