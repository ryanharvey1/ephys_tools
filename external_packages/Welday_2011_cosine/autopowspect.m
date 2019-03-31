function powspect = autopowspect(maxfreq, numfreqs, blocklen, autohisto, speedbins, balcounts)

%% maxfreq -- compute power spectrum from 0 - maxfreq
%% numfreqs -- number of frequency bins in the spectrum
%% blocklen -- length of each data block (in seconds)
%% ratehisto -- rate histogram data in matrix format (each row is the histogram for a movement episode of length blocklen)
%% speedbins -- an ordered list (length is same as number of rows in ratehisto) of which bin of the speed distribution each movement episide falls into
%% balcounts -- the counts for each speed bin in the balanced speed distribution (these define the "speed quotas" for the spectral analysis)

%% returns the power spectrum in 'powspect'

global iterat
f = (numfreqs/blocklen)*(0:(2^18))/2^19;  %frequency bins of the power spectrum
peakwidth = 1.5; %bandwidth (in Hz) on either side of the theta peak across which to intgrate for computing expected frequency value

nonzerobins=find(balcounts>0); %find the speed bin numbers with non-zero quotas

powspect=zeros(1,2^19); %initialize the power spectrum to zero for cumulative averaging
    
for i=1:100 %we will run the FFT a total of 100 times, randomly refilling the speed quotas each time
    
    ssignal=[]; %variable for accumulating the speed-balanced rate histogram on each pass through the data
    for j=nonzerobins %loop through the bins of the speed distribution which have nonzero quotas
        %%%first, try to fill the speed quota by sampling rows from the rate histogram WITHOUT replacement
        bsignal=[]; %variable for accumulating the rows of the rate histogram that fall into speed bin j
        thisbinepochs=autohisto(find(speedbins==j),:); %extract the rows of the autocorrelogram with running speeds that fall into speed bin j
        bsignal=[bsignal; thisbinepochs]; %accumulate the rows
        %%%second, if necessary, fill what remains of the speed quota by resampling rows from the rate histogram as little replacement as possible
        while (size(bsignal,1)<balcounts(j)) %while the quota for speed bin j has not yet been filled
            randthisbin=[thisbinepochs rand(size(thisbinepochs,1),1)]; %append a column of random numbers to the accumulated rate histogram
            randthisbin=sortrows(randthisbin,size(randthisbin,2)); %use the random number column to randomly order the rows
            randthisbin=randthisbin(:,1:size(thisbinepochs,2)); %remove he random number column now that its job is done
            if (size(bsignal,1)+size(randthisbin,1))<=balcounts(j) %if we have not yet exceeded our desired speed quota
                bsignal=[bsignal; randthisbin]; %then add all of the randomly ordered autocorrelogram rows to our accumulating total
            else %otherwise
                bsignal=[bsignal; randthisbin(1:(balcounts(j)-size(bsignal,1)),:)]; %add only as many rows as we need to fill the quota and be done
            end
        end
        ssignal=[ssignal; bsignal]; %accumulate the rows from speed bin j into the balanced rate histogram
    end
    
    Y=fft(sum(ssignal),2^19);
    Pyy = Y.* conj(Y) / 2^19;      
    powspect = powspect+Pyy/100;
    
 end





