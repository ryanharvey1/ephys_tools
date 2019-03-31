% walk
close all;
x=100;
y=100;

xmin=-500;
xmax=500;
ymin=-500;
ymax=500;

Xedges=xmin:30:xmax;
Yedges=ymin:30:ymax;

newx=x;
newy=y;
xout=x;
yout=y;
while true
    while true
        newx=randi([xout(end)-randi([0 50]),xout(end)+randi([0 50])]);
        if newx>xmin && newx<xmax
            break
        end
    end
    xout=[xout;newx];
    
    while true
        newy=randi([yout(end)-randi([0 50]),yout(end)+randi([0 50])]);
        if newy>ymin && newy<ymax
            break
        end
    end
    yout=[yout;newy];
    
    [N,~,~] = histcounts2(xout,yout,Xedges,Yedges);
    RC=reshape(N,[],1)';
    [~,nCells]=size(RC);
    if sum(RC,2).^2./(nCells*sum(RC.^2,2))>.80
        break
    end
end
x=smooth(xout,25);
y=smooth(yout,25);
% figure;comet(x,y)
[occmat,~,~] = histcounts2(y,x,Yedges,Xedges);
% figure;imagesc(occmat);colormap jet; axis off; box off;axis xy

RBS=false(1,length(x));
RBS(1:500)=true; 
spk=RBS(randperm(numel(RBS)))';

figure;
plot(x,y,'k');hold on
scatter(x(spk),y(spk),'filled','r')

[spkmat,~,~] = histcounts2(y(spk),x(spk),Yedges,Xedges);
ratemap=spkmat./occmat;

filtWidth = [5 5]; filtSigma = 1;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
ratemap = nanconv(ratemap,imageFilter, 'nanout');

figure;
imAlpha=ones(size(ratemap));
imAlpha(isnan(ratemap))=0;
imagesc(ratemap,'AlphaData',imAlpha);
colormap jet; axis off; box off;axis xy



frames=[linspace(0,length(x)/30,length(x))',x,y,spk];
spkts=frames(spk,1);

xmin=-500;
xmax=500;
ymin=-500;
ymax=500;


% repmat(xmin,1001,1)
% repmat(xmax,1001,1)
% repmat(ymin,1001,1);

cues=[[[xmin:xmax]',repmat(ymax,1001,1)];[repmat(xmax,1001,1),[ymax:-1:ymin]'];[[xmax:-1:xmin]',repmat(ymin,1001,1)];[repmat(xmin,1001,1),[ymin:ymax]']];

cues=cues(1:10:end,:);
% figure;plot(cues(:,1),cues(:,2),'.k')

clear p r
for i=1:length(cues)
    [p(i,1),r(i,1)]=ECD(frames,[],spkts,cues(i,:));
end

