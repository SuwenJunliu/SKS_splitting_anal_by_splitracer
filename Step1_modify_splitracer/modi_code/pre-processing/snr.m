function [snr]  =  snr(trace,dt,tw_noise,tw_signal)

% determine signal-to-noise ration
% usage:
% trace : time series, 100 s
% dt: sample rate in time steps
% tw_noise: length of noise window
% tw_signal: length of signal window

% Copyright 2016 M.Reiss and G.Rümpker

%modi by junliu
% sometimes the signal will show at few second before sks phase time

delay = 3;


% assume that signal starts at center of trace!
N = length(trace);
N_half = int32(N/2);

N_half = N_half - int32(delay/dt);

% determine beginning of noise window from tw_noise
N_noise = N_half-int32(tw_noise/dt) ;
% determine end of signal window from tw_signal
N_signal = N_half+int32(tw_signal/dt) ;


% max amp?
% calculate noise value
a_noise = 0;
for n = N_noise:N_half
    if a_noise < abs(trace(n))
        a_noise = abs(trace(n));
    end
    %a_noise = a_noise+abs(trace(n));
end

% calculate signal value
a_signal = 0;
for n = N_half:N_signal
    if a_signal < abs(trace(n))
        a_signal = abs(trace(n));
    end    
    %a_signal = a_signal+abs(trace(n));
end

nl_noise = (N_half-N_noise+1);
nl_signal = (N_signal-N_half+1);

%a_noise = a_noise/double(nl_noise);
%a_signal = a_signal/double(nl_signal);

snr = a_signal/a_noise;

end

