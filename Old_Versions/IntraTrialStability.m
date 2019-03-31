function r=IntraTrialStability(mat,linear_track,track_length)

if isempty(mat)
    r=NaN;
    return
end

if isequal(linear_track,'yes')
    nBinsy=1;
    filtWidth=[1,5];
elseif isequal(linear_track,'no')
    nBinsy=round(track_length/3);
    filtWidth = [5 5];
end


split=round(size(mat,1)/2);
S{1}=mat(1:split,:);
S{2}=mat(split+1:end,:);
MinY=min(mat(:,3));MaxY=max(mat(:,3));MinX=min(mat(:,2));MaxX=max(mat(:,2));
nBinsx = round(track_length/2);
edges{1}=linspace(MinY,MaxY,nBinsy+1);
edges{2}=linspace(MinX,MaxX,nBinsx+1);

for i=1:2
    mat=S{i};
    Omatrix = hist3([mat(mat(:,6)==0,3),mat(mat(:,6)==0,2)],'Edges',edges);
    Omatrix(end,:) = [];
    Omatrix(:,end) = [];
    occ = Omatrix/30;
    Smatrix = hist3([mat(mat(:,6)==1,3),mat(mat(:,6)==1,2)],'Edges',edges);
    Smatrix(end,:)=[];
    Smatrix(:,end)=[];
    FilledRateMatrix=Smatrix./occ;
    FilledRateMatrix(isinf(FilledRateMatrix))=0;
    imageFilter=fspecial('gaussian',filtWidth,1);
    map{i}=nanconv(FilledRateMatrix,imageFilter, 'nanout');
end
map1=map{1};
map2=map{2};
map1(isnan(map1))=0;
map2(isnan(map2))=0;
r=corrcoef(map1,map2);
r=r(2,1);
end