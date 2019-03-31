
% klustakwik_cleanup
files=[dir( '**/*.clu.1');dir( '**/*.fet.1');dir( '**/*.klg.1');dir( '**/*.model.1')];
disp([num2str(sum([files.bytes])/1000000000),'GB'])
toremove=strcat({files.folder}',filesep,{files.name}');
for i=1:size(files,1)
    delete(toremove{i});
end
