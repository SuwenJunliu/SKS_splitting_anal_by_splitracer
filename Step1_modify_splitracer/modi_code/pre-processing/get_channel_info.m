function [chan_info] = get_channel_info(txtfile,station)
% function to read infos about station, i.e. channel names & misalignment
% c. Frederik Link & Miriam Christina Reiss, 2019

% read channel file
fidCI = fopen(txtfile,'r');
formatSpec = '%s';
N = 17;
C_header = textscan(fidCI,formatSpec,N,'Delimiter','|');
C = textscan(fidCI,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',' ');
fclose(fidCI);
channel = C{1,4};
azimuth_txt = C{1,9};
dip_txt = C{1,10};
startdates = C{1,16};
enddates = C{1,17};

% find all different times spans in which the station settings (i.e. change
% in sensor / data logger / re alignment were changed
for j = 1:length(startdates)
    azimuth(j) = str2double(char(azimuth_txt{j}));
    dip(j) = str2double(char(dip_txt{j}));
    startvec{j} = datevec(char(startdates{j}(1:10)));
    startnum(j) = datenum(startvec{j});
    if isempty(enddates{j})
        endvec{j} = datevec(now);
    else
        endvec{j} = datevec(char(enddates{j}(1:10)));
    end
end

% sort those dates & assign channels

% get Z components
Z = contains(channel,'Z');
[row,~] = find(Z);

[fin_dates, ~, ~] = unique(startdates(row));

% check for all three channels,
jj = 0;
for j = 1:length(fin_dates)
    
    % get other channels
    check_date = strcmp(startdates,fin_dates(j));
    if sum(check_date)==3
        dp = dip(check_date);
        az = azimuth(check_date);
        ch = channel(check_date);
        % check if there is previous information
    elseif exist('chan_info','var') 
        % if there are one Z channel and one other
        ch2 = channel(check_date); 
        if length(ch2) == 2
            % check which other channel exists and only use previous
            % information for missing channel, takes previous ch, az, dp
            if any(contains(ch2,'N')) || any(contains(ch2,'1')) 
                dp_tmp = dip(check_date);
                az_tmp = azimuth(check_date);
                dp(1) = dp_tmp(1);
                az(1) = az_tmp(1);
                ch(1) = ch2(1);
                dp(3) = dp_tmp(2);
                az(3) = az_tmp(2);
                ch(3) = ch2(2);
            elseif any(contains(ch2,'E')) || any(contains(ch2,'2'))
                dp_tmp = dip(check_date);
                az_tmp = azimuth(check_date);
                dp(2) = dp_tmp(1);
                az(2) = az_tmp(1);
                ch(2) = ch2(1);
                dp(3) = dp_tmp(2);
                az(3) = az_tmp(2);
                ch(3) = ch2(2);
            end
            % if there is only a Z component, take previous settings from ch,az,dp
        else 
            dp(3) = dip(check_date);
            az(3) = dip(check_date);
            ch(3) = ch2;
        end
    else 
        disp('no channel information for station / dates found')
        continue
    end
    
    %    checks whether Z is vertical,
    indx = 1:length(dp);
    [~,I] = max(abs(dp));
    dipZ = dp(I);
    if dipZ ~= -90
        continue
    end
    
    % check that N
    % &E components are aligned 90? to each other and check whether they ares
    % switched
    channel_Z = ch(I);
    az_oi = az(indx~=I);
    ch_oi = ch(indx~=I);
    az_flag = round(mod(az_oi(1)-az_oi(2),360));
    switch az_flag
        case 90
            Nidx = 2;
            Eidx = 1;
            
        case 270
            Nidx = 1;
            Eidx = 2;
        otherwise
            startv = startvec(check_date); sv = startv{end};
            endv = endvec(check_date); ev = endv{end};
            date_str = [char(datetime(sv)),' - ',char(datetime(ev))];

            w = warndlg(['Horizontal components are not perpendicular for',...
                'the time period from ',date_str, ...
                'Please check channel info file if you want to make ',...
                'any manuel changes'],['Warning for station: ',char(station)]);
            drawnow
            waitfor(w);
            continue
    end
    azN = az(Nidx);
    azE = az(Eidx);
    channel_N = ch_oi(Nidx);
    channel_E = ch_oi(Eidx);
    
    % calculate whether components are misaligned
    cor_deg = azN;
    if cor_deg > 180
        cor_deg = cor_deg-360;
    elseif cor_deg < -180
        cor_deg = cor_deg+360;
    end
    if cor_deg == 0
        rot_flag = 0;
    else
        rot_flag = 1;
    end
    if length(channel_Z) > 1
        channel_Z = char(channel_Z{1});
    end
    if length(channel_N) > 1
        channel_N = char(channel_N{1});
    end
    if length(channel_E) > 1
        channel_E = char(channel_E{1});
    end
    
    jj = jj +1;
    
    startv = startvec(check_date);
    endv = endvec(check_date);
    
    % save infos
    chan_info(jj).start_vec = startv{end}; % take time vec from Z
    chan_info(jj).end_vec = endv{end};
    chan_info(jj).channel_Z = channel_Z;
    chan_info(jj).channel_N = channel_N;
    chan_info(jj).channel_E = channel_E;
    chan_info(jj).rot_flag = rot_flag;
    chan_info(jj).cor_deg = cor_deg;
end


end