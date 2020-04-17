% 
% x=zeros(70,30)
% x(round(xcoords),round(ycoords))=chanMap'
% 
% for i=1:32
%     x(round(ycoords(i)),round(xcoords(i)))=chanMap(i)
% 
% end

figure
map=zeros(200);
imagesc(map);
hold on
grid on
set(gca,'GridColor',[1 1 1])
xlabel('probe x [\mum]','interpreter','Tex')
ylabel('probe y [\mum]','interpreter','Tex')

i=1;
j=1;
while i==1
    [xcoords(j),ycoords(j)]=ginput(1);
    xcoords=round(xcoords,-1);
    ycoords=round(ycoords,-1);
    map(ycoords(j),xcoords(j))=1;
    imagesc(map);
    j=j+1;
end



