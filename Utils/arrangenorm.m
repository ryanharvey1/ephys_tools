function [matout]=arrangenorm(mat)
% FOR DIAG POP VECTOR
[~,I]=max(mat,[],2);
[~,I2]=sort(I);
matout=mat(I2,:);
for i=1:size(matout,1)
    matout(i,:)=rescale(matout(i,:),0,1);
end
end