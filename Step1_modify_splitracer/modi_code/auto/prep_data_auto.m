function prep_data_auto(sel_data)

% initial preparation of data before further analysis

% Copyright 2021 F.Link, M.Reiss and G.Rümpker 


% get selected stations
[selected_stations,stations] = read_stationfile(sel_data);

% create new directory
dir_save = [sel_data.set_folder,'/pre-processing/'];
if ~exist(dir_save,'dir')
    mkdir(dir_save)
end


st_ind = strcmp(char(selected_stations),stations.name);
disp(char(strcat({'Processing station '},selected_stations,{' '},char(stations.nc(st_ind)))))

% check if pre-processed file already exists
file_save = [dir_save,char(selected_stations),'.mat'];
last_event = [dir_save,char(selected_stations),'_last_event.mat'];

if exist(file_save,'file') && exist(last_event,'file')
    disp('A pre-processed .mat-file for this station already exists. Skipping station.')
    return
end

% get channel information
[chan_info] = get_channel_info_auto(strcat(sel_data.work_dir, ...
    '/mseed_files/', char(selected_stations), '_', ...
    char(stations.nc(st_ind)),'/channel_info.txt'),selected_stations);

%read log file containing successfully downloaded data
fileID2 = fopen([sel_data.work_dir,'/mseed_files/data_request_',...
    char(selected_stations),'_',char(stations.nc(st_ind))]);
C2 = textscan(...
    fileID2,'%f %f %f %f %f %f %f %f %f %f %f %f %s %s');
fclose(fileID2);
files = C2{1,14};

prev_n = 1;
ec = 0;
% loop around all files per station
for n =  prev_n:length(files)
    
    % read event from text file & compare with channel info to find out
    % whether miniseed must be rotated because of misalignment
    [event] = read_event_file(C2,n);
    eventdatevec = datevec(...
        [num2str(event.year) '-' num2str(event.month) '-' num2str(event.day)]);
    eventnum = datenum(eventdatevec);
    if length(chan_info) > 1
        found_ch = 0;
        for j = 2:length(chan_info)
            if eventnum<datenum(chan_info(j).start_vec)
                choi = j-1;
                found_ch = 1;
                break
            end
        end
        if ~found_ch
            choi = length(chan_info);
        end
    else
        choi = 1;
    end

    % read miniseed
    disp(char(files(n)));
    file_read = strcat(char(files(n)));
    file_mseed = extract_adv_auto_sac(file_read, chan_info(choi));

    % check if more than 2 components exist
    if length(file_mseed(1,:))>2

    if ~isempty(file_mseed(1).comp) && ~isempty(file_mseed(2).comp) && ...
            ~isempty(file_mseed(3).comp) 
        % check if components consist of more than 30 minutes
        no_samples = 1800 / file_mseed(1).dt;

        if sum(size(file_mseed(1).amp))>no_samples && ...
                sum(size(file_mseed(2).amp))>no_samples && ...
                sum(size(file_mseed(3).amp))>no_samples

            % resample data
            [event] = re_sample(event,file_mseed,0.05);

            % rotate if necessary
            if chan_info(choi).channel(file_mseed(1).chID).rot_flag
                [event.n_amp,event.e_amp] = rot_az(event.n_amp,...
                    event.e_amp, chan_info(choi).channel(file_mseed(1).chID).cor_deg);
            end

            %get travel times
            [event.phases] = get_tt(event.dist,event.depth);

            if isstruct(event.phases) == 0
                disp('travel times could not be calculated')
                continue
            end

            % calculate absolute travel times
            for iF = 1:length(event.phases)
                event.phases(iF).tt_abs = ...
                    seconds(event.phases(iF).tt) + event.origin_time;
            end

            %% first quality check
            event = rc_phases(event,sel_data);
            to_save = 0;

            % check if phases fulfill set criteria
            for i_phase = 1:length(event.phases)
                if event.phases(i_phase).snr > sel_data.snr
                    if isfield(event.phases(i_phase),'cordeg')
                    to_save = to_save +1;
                    event.phases_to_analyze(to_save) = i_phase; 
                    end
                end
            end
            if to_save == 0
                disp('SNR too low')
            end

            % save phase
            if to_save > 0
                ec = ec + 1;
                new_name = ['event' int2str(ec)];
                data.(new_name) = event;
            end

            clear event
        else
            disp(['individual components of this event',...
                ' are too short for analysis'])
        end
    else
        disp('event does not have three components')
    end
    else
        disp('event does not have three components')
    end
end
if exist('data','var') == 1
    save([dir_save,char(selected_stations),'.mat'],'data');
    save([dir_save,char(selected_stations),'_last_event.mat'],'n');
else
    disp('there are no events for this station')
end

clear data n


end