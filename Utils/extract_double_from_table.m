function [mat,varnames]=extract_double_from_table(tab)
fields= fieldnames(tab);
    mat=[];
    varnames=[];
    for i=1:size(tab,2)
        if isa([tab.(fields{i}){:}], 'double') || isa([tab.(fields{i}){:}], 'logical')
            mat=[mat,cell2mat(tab.(fields{i}))];
            varnames=[varnames,fields(i)];
        end
    end
end