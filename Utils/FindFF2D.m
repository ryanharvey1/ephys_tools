function [field,fieldWidth]=FindFF2D(Ratemap)
% FindFF2D Finds firing field in a firing map, i.e. the connex area around the peak,
% where the firing rates are above a given threshold (20% of max bin by default).
% Fields will have at least 10 bins
%
% Input:
%              Ratemap:     Binned matrix of rates
%
% Output:
%                field:     Logical index of field location
%           fieldWidth:     Width of field in bins
%
% Ryan E. Harvey & Laura B., July 18th 2017

% take care of NaN maps
if sum(sum(isnan(Ratemap)))==(size(Ratemap,1)*size(Ratemap,2))
    field=zeros(length(Ratemap));
    fieldWidth=0;
    return
end

if nansum(nansum(Ratemap))==0
    fieldWidth=0;
    field=ones(length(Ratemap));
    return
end

% pad maps
Ratemap=padarray(Ratemap,[2 2]);

% sort rates
r = sort(reshape(Ratemap,[],1),'descend'); r(isnan(r))=[];
for ii=1:length(r)
    M=r(ii); % recursive call if field is too small
    if M==0 % find really small fields
        if exist('field','var')
        field=field(3:end-2,3:end-2); 
        field2=field; field2(~any(field2,2),:)=[]; field2(:,~any(field2,1))=[];
        fieldWidth=max([sum(diag(field2,0)),size(field2,1),size(field2,2)]);
        else
            field=NaN; fieldWidth=NaN;
        end
        return 
    end
    % FIND THRESHOLD AND MAX
    Thres=M*0.20; [row,col]=find(Ratemap==M); row=row(1); col=col(1);
    
    % CREATE EMPTY LOGICAL MAP
    field=zeros(size(Ratemap)); field(row,col)=1;
    
    % CHECK INITAL SURROUNDINGS WITH 3X3 BOX
    boxx=Ratemap(row-1:row+1,col-1:col+1); 
    result=boxx>Thres;
    field(row-1:row+1,col-1:col+1)=result;
    
    % CHECK ALL SURROUNDINGS FOR CONTIGUOUS BINS OVER THRESHOLD
    lagField=[];
    noGrowth=1;
    [row,col]=find(field==1);
    for j=1:length(Ratemap)
        if sum(sum(lagField))==sum(sum(field))
            noGrowth=noGrowth+1;
            if noGrowth==3
                field2=field; field2(~any(field2,2),:)=[]; field2(:,~any(field2,1))=[];
                fieldWidth=max([sum(diag(field2,0)),size(field2,1),size(field2,2)]);
                if sum(sum(field))>=10 || sum(sum(field))==sum(sum(Ratemap>0)) % have at least 10 bins in field
                    field=field(3:end-2,3:end-2);
                    return
                else
                    break % if field is under 10 bins, redefine field using recursively
                end
            end
        end
        lagField=field;
        [row,col]=find(field==1);
        for i=1:length(row)
            boxx=Ratemap(row(i)-1:row(i)+1,col(i)-1:col(i)+1);
            result=boxx>Thres;
            lagField=field;
            field(row(i)-1:row(i)+1,col(i)-1:col(i)+1)=result;
        end
    end
end
end


