function [speed_dist, speedvec, ratehisto, autohisto] = block_speeds(M_d,Position_Speed,spikedata,pixpercm,blocksize,binsize,speed_edges, blocksec, trackhz)

%%to prepare for speed balancing of movement intervals, this function 
%%divides speed intervals into blocks of size 'blocksize'
%%and computes the running speed in each block

%%M_d is the vector of movement intervals in direction being analyzed(e_fints, ne_fints, etc.)
%%Position_Speed(:,1) contains sample speeds in pix per sample
%%Position_Speed(:,2) contains sample time stamps in sec
%%pixpercm is tracker resolution in pixels per cm 
%%blocksize is the number of position samples per blocksec
%%binsize is width of each time bin for the firing rate histogram
%%speed edges are the bin boundaries (in cm/s) for the running speed distribution
%%blocksec is the length (in sec) of each movement block that is being analyzed
%%trackhz is the sample rate of the position tracker in hz

%speed_dist returns the distribution of running speeds (binned by 'speed_edges') for the set 
%            of movement intervals of length 'blocksec' in the direction being analyzed
%speedvec returns 5 columns of data, with each row corresponding to a single movement episode of length 'blocksec':
%    col1: the mean speed of the movement block (in cm/s)
%    col2: start time of the block
%    col3: end time of the block
%    col4: row number of 'M_d' from which the interval was extracted
%    col5: the bin number of 'speed_dist' into which this movement episode was placed
%ratehisto returns matrix whose rows are the rate histogram of the theta cell during a single movement block
%autohisto returns matrix whose rows are the autocorrelogram of the theta cell during a single movement block

cmperpix=1/pixpercm;
halfsample=1/(trackhz*2);

speedvec=[];
ratehisto=[];
autohisto=[];
intvec=[];

tempints=M_d; 

for i=1:size(tempints,1) %%loop through the movement intervals
    idex=find((Position_Speed(:,2)>=tempints(i,1)) & (Position_Speed(:,2)<tempints(i,2))); %find indices into Position_Speed for the current interval
    if length(idex)>0  %if there is speed data for this interval in the time range being analyzed
        ispeeds=Position_Speed(idex,1); %get speeds at each position sample in the interval
        itimes=Position_Speed(idex,2); %get time stamps of each position sample in the interval
        leftover=mod(length(ispeeds),blocksize); %we must truncate this many position samples rom the interval so its length will be an integer multiple of 'blocksec'
        numblocks=(length(ispeeds)-leftover)/blocksize; %%integer number of blocks in this interval
        if (numblocks>=1) %if there is at least one block in this interval
            chopbeg=mean(ispeeds((1+leftover):length(ispeeds))*cmperpix*trackhz); %mean interval speed that will result if excess samples are removed from the beginning of the interval
            chopend=mean(ispeeds(1:(length(ispeeds)-leftover))*cmperpix*trackhz); %mean interval speed that will result if excess samples are removed from the end of the interval
            %%trim beg or end, whichever yields highest speed mean for the entire interval:
            if (chopbeg>=chopend) %trim samples from beginning of the interval if it yields the highest mean interval speed
                tempints(i,1)=itimes(1+leftover); %adjust the timestamp for the start of the interval
                for j=1:numblocks %add a row to 'speedvec' for each movement block in the interval 
                    newline=[mean(ispeeds((1+leftover+(j-1)*blocksize):(blocksize+leftover+(j-1)*blocksize))*cmperpix*trackhz) itimes(1+leftover+(j-1)*blocksize) itimes(blocksize+leftover+(j-1)*blocksize) i];
                    if ((newline(3)+halfsample)-(newline(2)-halfsample))>blocksec %if the block has enough position samples, add it to 'speedvec'
                        speedvec=[speedvec; newline];
                        ratechunk=histc(spikedata,[(newline(2)-halfsample):binsize:(newline(3)+halfsample)]); %also add a row to the 'ratehisto' matrix
                        ratehisto=[ratehisto; ratechunk(1:(length(ratechunk)-1))'];
                        spikedex=find(spikedata<(newline(3)+halfsample));  %also add a row to the 'autohisto' matrix
                        spikedex=find(spikedata(spikedex)>=(newline(2)-halfsample));
                        autochunk=autocorrelogram(spikedata(spikedex), (blocksec/256), blocksec);
                        autohisto=[autohisto; autochunk];
                    end
                end
            else %trim samples from end of interval if this yields the highest mean interval speed
                tempints(i,2)=itimes(length(itimes)-leftover);  %adjust the timestamp for the end of the interval
                for j=1:numblocks %add a row to 'speedvec' for each movement block in the interval 
                    newline=[mean(ispeeds((length(ispeeds)-leftover-(j-1)*blocksize-(blocksize-1)):(length(ispeeds)-leftover-(j-1)*blocksize))*cmperpix*trackhz) itimes(length(ispeeds)-leftover-(j-1)*blocksize-(blocksize-1)) itimes(length(ispeeds)-leftover-(j-1)*blocksize) i];
                    if ((newline(3)+halfsample)-(newline(2)-halfsample))>blocksec  %if the block has enough position samples, add it to 'speedvec'
                        speedvec=[speedvec; newline];
                        ratechunk=histc(spikedata,[(newline(2)-halfsample):binsize:(newline(3)+halfsample)]);  %also add a row to the 'ratehisto' matrix
                        ratehisto=[ratehisto; ratechunk(1:(length(ratechunk)-1))'];
                        spikedex=find(spikedata<(newline(3)+halfsample));  %also add a row to the 'autohisto' matrix
                        spikedex=find(spikedata(spikedex)>=(newline(2)-halfsample));
                        autochunk=autocorrelogram(spikedata(spikedex), (blocksec/256), blocksec);
                        autohisto=[autohisto; autochunk];
                    end
                end
            end
        end
    end
end

ratehisto=[ratehisto speedvec(:,2)]; %append a column of movement block starttimes to the right end of the rate histogram
ratehisto=sortrows(ratehisto,17); %make sure the rows of the rate histo are in correct temporal order, just in case
ratehisto=ratehisto(:,1:16); %remove the column of movement block starttimes from the rate histogram
autohisto=[autohisto speedvec(:,2)]; %append a column of movement block starttimes to the right end of the autocorrelograms
autohisto=sortrows(autohisto,514); %make sure the rows of autocorrelograms are in correct temporal order, just in case
autohisto=autohisto(:,1:513); %remove the column of movement block starttimes from the autocorrelogram
speedvec=sortrows(speedvec,2); %make sure that 'speedvec' is also in correct temporal order
speedvec(:,2)=speedvec(:,2)-halfsample; %adjust start time of speed interval to be in between position samples
speedvec(:,3)=speedvec(:,3)+halfsample; %adjust end time of speed interval to be in between position samples
[speed_dist,c]=histc(speedvec(:,1),speed_edges); %compute the distribution of running speeds for the direction being analyzed
speedvec=[speedvec c]; %append a column denoting which bin number of 'speed_dist' each movement block was placed into

    %%%%%detrrend the rate histogram
    numrows=size(ratehisto,1); %# of rows in the rate histogram
    dt=ratehisto'; %transpose the rate histogram
    dt=dt(:);      %convert to a vector
    bm=polyfit(0:(length(dt)-1),dt',1); %perform linear regression on the rate histogram across the entire session
    regline=0:bm(1):(bm(1)*(length(dt)-1)); %compute the regression line
    dt=dt-regline';  %detrend by substracting the regression line from the rate histogram
    dt=reshape(dt,16,numrows); %restore the rate histogram to matrix format
    ratehisto=dt'; %de-transpose

    
 function acor = autocorrelogram(spiketimes, binsize, timewidth)
   if ~isempty(spiketimes)
       
     if size(spiketimes,1)<size(spiketimes,2)
         spiketimes=spiketimes';
     end
         
     s2=[spiketimes; spiketimes];
     ISIs=[];
     
     for i=1:length(spiketimes)
         ISIs = [ISIs; spiketimes(i)-s2((i+1):(i+length(spiketimes)))];
     end
     
     ISIs=ISIs(find(~(ISIs==0)));
     if ~isempty(ISIs)
         acor = histc(ISIs,-timewidth:binsize:timewidth);
     else
         acor=[-timewidth:binsize:timewidth]*0;
     end
    
   else
       
     acor=[-timewidth:binsize:timewidth]*0;
       
   end
   
     if size(acor,1)>size(acor,2)
         acor=acor';
     end
    






