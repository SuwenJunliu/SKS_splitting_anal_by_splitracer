function [trace]= extract_adv_auto(file,chan_info)

% this function reads the mutilple sac file, double entries
% a low pass of 1 s is administered to avoid aliasing effects
% that's only substitude the module for reading mseed


% usage:
% file: mseed file

% Copyright 2022 Junliu Suwen @ IGGCAS



try
s = readsac([file '*']);
catch 
    disp('could not read mseed')
    trace = 0;
    return
end

% Check if correct station is in the header
for i = 1:length(s)
    if  ~strcmp(chan_info.net,s(i).KNETWK) || ~strcmp(chan_info.stat,s(i).KSTNM)
        disp('wrong file downloaded')
        trace = 0;
        return
    end
end



for i = 1:length(chan_info.channel)
    locID(i) = str2double(chan_info.channel(i).location);
    locID(isnan(locID)) = -1;
    len(i) = chan_info.channel(i).chnum;
end
[locID_u,J] = unique(locID);
for i = 1:length(locID_u)
    loc_flag(i) = len(i)==3;
end
if isempty(locID_u(loc_flag))
    disp('Less than three components measured in this period')
    trace = 0;
    return
end
[~,K] = min(locID_u);
locstr = chan_info.channel(J(K)).location;
chanstr = [chan_info.net ':' chan_info.stat ':' char(locstr) ':'];

foundflag = 0;
for i_comp=1:length(s)
    
    clear comp_name
    clear comp
    clear trace0
    clear t
    
    comp_name = s(i_comp).KCMPNM;

    
    trace0 = s(i_comp).DATA1;
    trace0 = trace0-mean(trace0);
    dt = s.DELTA;
    %trace0 = buttern_low(trace0',6,1,dt);
    trace0 = buttern_filter(trace0',6,0.01,1,dt);
    [year,month,day] = cal_day(s(i_comp).NZYEAR,s(i_comp).NZJDAY);
    starttime = datenum(year,month,day,s(i_comp).NZHOUR,s(i_comp).NZMIN,s(i_comp).NZSEC + s(i_comp).B);
    endtime =  datenum(year,month,day,s(i_comp).NZHOUR,s(i_comp).NZMIN,s(i_comp).NZSEC + s(i_comp).E);
    t = linspace(starttime,endtime,s(i_comp).NPTS);


    % use appropriate channel name
    comp_name = char(comp_name);
    
    if contains(comp_name,"L")
        found_flag = 1;
    elseif  contains(comp_name,"B")
        found_flag = 0;
    end
    
    
    if ~found_flag
    for ii = 1:length(chan_info.channel)
        compref_Z = [chanstr char(chan_info.channel(ii).channel_Z)];
        compref_N = [chanstr char(chan_info.channel(ii).channel_N)];
        compref_E = [chanstr char(chan_info.channel(ii).channel_E)];
    if contains(compref_Z,comp_name)
        trace(3).comp = 'Z';
        trace(3).amp = trace0;
        trace(3).time = t;
        trace(3).dt = dt;
        trace(3).chID = ii;
        found_flag = 1;
        L = ii;
    elseif contains(compref_E,comp_name)
        trace(1).comp = 'East';
        trace(1).amp = trace0;
        trace(1).time = t;
        trace(1).dt = dt;
        trace(1).chID = ii;
        found_flag = 1;
        L = ii;
    elseif contains(compref_N,comp_name)
        trace(2).comp = 'North';
        trace(2).amp = trace0;
        trace(2).time = t;
        trace(2).dt = dt;
        trace(2).chID = ii;
        found_flag = 1;
        L = ii;
    end
    end
    else
        for ii = 1:length(chan_info.channel)
        compref_Z = [chanstr char(chan_info.channel(ii).channel_Z)];
        compref_N = [chanstr char(chan_info.channel(ii).channel_N)];
        compref_E = [chanstr char(chan_info.channel(ii).channel_E)];
        if contains(compref_Z,comp_name)
            trace(3).comp = 'Z';
            trace(3).amp = trace0;
            trace(3).time = t;
            trace(3).dt = dt;
            trace(3).chID = ii;
        elseif contains(compref_E,comp_name)
            %keyboard
            trace(1).comp = 'East';
            trace(1).amp = trace0;
            trace(1).time = t;
            trace(1).dt = dt;
            trace(1).chID = ii;
        elseif contains(compref_N,comp_name)
            trace(2).comp = 'North';
            trace(2).amp = trace0;
            trace(2).time = t;
            trace(2).dt = dt;
            trace(2).chID = ii;
        end
        end
    end
end


if ~exist('trace','var')
    trace(1).comp = 'NotReadable';
    trace(1).amp = 0;
    trace(1).time = 0;
    trace(1).dt = 0;
    trace(1).chID = 0;
    return
end


end