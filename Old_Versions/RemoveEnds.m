function [ matrix1,matrix2,track_length ] = RemoveEnds( pixelDist,matrix1,matrix2,track_length,CM2Remove )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Pix2remove=how many pixels in 5cms of track
pix2remove=CM2Remove/pixelDist; 
% FIND MAX AND MIN IN X DIM
x_min=min(matrix1(:,2)); x_max=max(matrix1(:,2));
% FIND VALUE OF PIXEL CUT OFF
lowercut=x_min+pix2remove; highercut=x_max-pix2remove;
% CUT OUT PIXELS
matrix1=matrix1((matrix1(:,2)>lowercut & matrix1(:,2)<highercut),:);
matrix2=matrix2((matrix2(:,2)>lowercut & matrix2(:,2)<highercut),:);
% ADJUST NEW TRACK LENGTH
track_length=track_length-CM2Remove*2;
end

