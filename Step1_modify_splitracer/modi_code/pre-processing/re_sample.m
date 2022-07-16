function [event]=re_sample(event,trace,new_sr)
% this function resamples event traces
% the interp1 function is used to have the data points starting exactly at
% origin time

% Copyright 2016 M.Reiss and G.Rümpker

% read origin time
date_val = [event.year,event.month,event.day,event.hour,event.min,event.sec];
event.origin_time = datetime(date_val);
event.ot_num=datenum(event.origin_time);

% resample
dt=datenum(seconds(new_sr));

% define new time vector
start_t = event.ot_num ; 
end_t =event.ot_num + datenum(seconds(2300)) ;
new_dt =start_t:dt:end_t;

for n=1:length(trace)
    trace(n).new_time = new_dt;
    trace(n).new_amp = interp1(trace(n).time,trace(n).amp,...
        trace(n).new_time,'linear','extrap');
end

event.sr = new_sr;
event.z_amp = trace(3).new_amp;
event.n_amp = trace(2).new_amp;
event.e_amp = trace(1).new_amp;
event.time = trace(n).new_time;

event = rmfield(event,'year');
event = rmfield(event,'month');
event = rmfield(event,'day');
event = rmfield(event,'hour');
event = rmfield(event,'min');
event = rmfield(event,'sec');

end