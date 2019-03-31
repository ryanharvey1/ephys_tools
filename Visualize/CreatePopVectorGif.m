function popvecgif_fig=CreatePopVectorGif(mapright_1C,mapleft_1C,mapright_2C,mapleft_2C)
% Create popvector gif
% attempt 2
pathtosave={'/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/PAE Project/Presentations/T32_JournalClub/PopVecGif_pae'};

tempmatTL=NaN(size(mapright_1C));
tempmatTR=NaN(size(mapleft_1C));
tempmatBL=NaN(size(mapright_2C));
tempmatBR=NaN(size(mapleft_2C));

popvecgif_fig=figure;
popvecgif_fig.Color=[1 1 1];
popvecgif_fig.OuterPosition=[1 1 1082 1058];


subplot(2,2,1)
imagesc(tempmatTL);colormap jet;box off;axis off
title('Right Run')
set(gca,'FontSize', 25)


subplot(2,2,2)
imagesc(tempmatTR);colormap jet;box off;axis off
title('Left Run')
set(gca,'FontSize', 25)


subplot(2,2,3)
imagesc(tempmatBL);colormap jet;box off;axis off


subplot(2,2,4)
imagesc(tempmatBR);colormap jet;box off;axis off

nframe=1;
print(popvecgif_fig,'-djpeg',[pathtosave{1},filesep,['popvecgif_fig',num2str(nframe),'.jpeg']]);nframe=nframe+1;



for i1=1:length(tempmatTL)
    tempmatTL(i1,:)=mapright_1C(i1,:);
    subplot(2,2,1)
    imagesc(tempmatTL);colormap jet;box off;axis off
    title('Right Run')
    set(gca,'FontSize', 25)
    print(popvecgif_fig,'-djpeg',[pathtosave{1},filesep,['popvecgif_fig',num2str(nframe),'.jpeg']]);nframe=nframe+1;
    
end
for i4=1:length(tempmatBR)
    tempmatBR(i4,:)=mapleft_2C(i4,:);
    subplot(2,2,4)
    imagesc(tempmatBR);colormap jet;box off;axis off
    print(popvecgif_fig,'-djpeg',[pathtosave{1},filesep,['popvecgif_fig',num2str(nframe),'.jpeg']]);nframe=nframe+1;
end
for i2=1:length(tempmatTR)
    tempmatTR(i2,:)=mapleft_1C(i2,:);
    subplot(2,2,2)
    imagesc(tempmatTR);colormap jet;box off;axis off
    title('Left Run')
    set(gca,'FontSize', 25)
    print(popvecgif_fig,'-djpeg',[pathtosave{1},filesep,['popvecgif_fig',num2str(nframe),'.jpeg']]);nframe=nframe+1;
end
for i3=1:length(tempmatBL)
    tempmatBL(i3,:)=mapright_2C(i3,:);
    subplot(2,2,3)
    imagesc(tempmatBL);colormap jet;box off;axis off
    print(popvecgif_fig,'-djpeg',[pathtosave{1},filesep,['popvecgif_fig',num2str(nframe),'.jpeg']]);nframe=nframe+1;
end

% REARANGE RATEMAPS IN DIAG
[mapright_1C,mapleft_1C]=arrange(mapright_1C,mapleft_1C);
[mapleft_2C,mapright_2C]=arrange(mapleft_2C,mapright_2C);

for i1=1:length(tempmatTL)
    tempmatTL(i1,:)=mapright_1C(i1,:);
    subplot(2,2,1)
    imagesc(tempmatTL);colormap jet;box off;axis off
    title('Right Run')
    set(gca,'FontSize', 25)
    print(popvecgif_fig,'-djpeg',[pathtosave{1},filesep,['popvecgif_fig',num2str(nframe),'.jpeg']]);nframe=nframe+1;
end
for i4=1:length(tempmatBR)
    tempmatBR(i4,:)=mapleft_2C(i4,:);
    subplot(2,2,4)
    imagesc(tempmatBR);colormap jet;box off;axis off
    print(popvecgif_fig,'-djpeg',[pathtosave{1},filesep,['popvecgif_fig',num2str(nframe),'.jpeg']]);nframe=nframe+1;
end
for i2=1:length(tempmatTR)
    tempmatTR(i2,:)=mapleft_1C(i2,:);
    subplot(2,2,2)
    imagesc(tempmatTR);colormap jet;box off;axis off
    title('Left Run')
    set(gca,'FontSize', 25)
    print(popvecgif_fig,'-djpeg',[pathtosave{1},filesep,['popvecgif_fig',num2str(nframe),'.jpeg']]);nframe=nframe+1;
end
for i3=1:length(tempmatBL)
    tempmatBL(i3,:)=mapright_2C(i3,:);
    subplot(2,2,3)
    imagesc(tempmatBL);colormap jet;box off;axis off
    print(popvecgif_fig,'-djpeg',[pathtosave{1},filesep,['popvecgif_fig',num2str(nframe),'.jpeg']]);nframe=nframe+1;
end

% ARRANGE BINS RIGHT
    function [mat1,mat2]=arrange(mat1,mat2)
        [~,I]=max(mat1,[],2);
        [~,I2]=sort(I);
        mat1=mat1(I2,:);
        mat2=mat2(I2,:);
    end

end

