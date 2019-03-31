function MyPointer(name)
fname= ['/u12/antsiro/matlab/draft/' name '.map'];
if exist(fname,'file')
    map = load(fname);
    set(gcf,'PointerShapeCData',map);
    set(gcf,'Pointer','custom');
else
    fprintf('cursor %s doesnot exist',name);
end