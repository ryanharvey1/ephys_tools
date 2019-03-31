 function mival = mi_measure(phase_sig, amp_sig)
% function mival = mi_measure(phase_sig, amp_sig)
%
% Returns a value for the MI measure calculated between two signals.
% (Functionality to deal with multiple trials will be added soon)
%
% INPUTS:
%
% phase_sig - the instantaneous phase values for a signal which has been
% filtered for a lower, modulating frequency, passed as a column vector
%
% amp_sig - the amplitude values for a signal which has been filtered for a
% higher, modulated frequency, passed as a column vector 
%
% Author: Angela Onslow, May 2010

num_trials = size(phase_sig, 2);

for count = 1:num_trials
    
    %Create composite signal
    z = amp_sig(:,count).*exp(1i*phase_sig(:,count));
    m_raw(count) = mean(z);  %Compute the mean length of composite signal.
        

    mival(count,1) = abs((m_raw(count)));
end

if num_trials > 1
    mival = mean(mival);
end

