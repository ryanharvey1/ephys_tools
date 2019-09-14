function [spikes, times]=load_pxd_data(fileroot, cell)
% load <fileroot>_pxdt.txt and <fileroot>_c<cell>_pxds.txt datafiles from single cell recording studies. 
% E.g.:
% [spikes, times]=load_pxd_data('pc', 1);
% [spikes, times]=load_pxd_data('hd', 1);
% [spikes, times]=load_pxd_data('tpd', 1);
% spikes and times arrays have dimensions n_dir_bins, n_x_bins, n_y_bins (n_x_bins = n_y_bins)

tfile=[fileroot '_pxdt.txt'];
sfile=[fileroot '_c' num2str(cell) '_pxds.txt'];

tfn=fopen(tfile,'r');
sfn=fopen(sfile,'r');

% read initial dim = n_pxn_pxn_d: 
s=fscanf(tfn,'%c',12);
tdim=fscanf(tfn,'%d',3);
s=fscanf(sfn,'%c',12);
sdim=fscanf(sfn,'%d',3);

if tdim~=sdim
   error('Dimensions do not match');
else
   dim=[0 0 0];
   dim(1)=tdim(1);
   dim(2)=tdim(2);
   dim(3)=tdim(3);
end

spikes=fscanf(sfn,'%d');
times=fscanf(tfn,'%f');

spikes=reshape(spikes,dim(3),dim(2),dim(1));
times=reshape(times,dim(3),dim(2),dim(1));

for i=1:dim(3)
    spikes(i,:,:) = squeeze(spikes(i,:,:))';
    times(i,:,:) = squeeze(times(i,:,:))';
end

fclose(sfn);
fclose(tfn);
