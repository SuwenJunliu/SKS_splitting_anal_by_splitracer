function prep_data(sel_data)

% initial preparation of data before further analysis

% Copyright 2016 M.Reiss and G.Rümpker
% altered 2019

% get selected stations
[selected_stations,stations] = read_stationfile(sel_data);

% create new directory
dir_save = [sel_data.set_folder,'/pre-processing/'];
mkdir(dir_save)

% loop of selected stations
for i = 1:length(selected_stations)
    
    st_ind = strcmp(char(selected_stations(i)),stations.name);
    
    % check if pre-processed file already exists
    file_save = [dir_save,char(selected_stations(i)),'.mat'];
    last_event = [dir_save,char(selected_stations(i)),'_last_event.mat'];
    
    if exist(file_save,'file') && exist(last_event,'file')
        outf = warning_3buttons(...
            'A pre-processed .mat-file for this station already exists',...
            'append','overwrite',char(selected_stations(i)));
        if outf == 1
            % load previous last event index
            load(last_event); prev_n = n+1;  clear('n');
            % load saved phases
            load([dir_save,char(selected_stations(i)),'.mat'])
            ec = length(fieldnames(data));  
        elseif    outf == 2
            prev_n = 1;
            ec = 0;
        else
            continue
        end
    else 
        prev_n = 1;
        ec = 0;
    end
    
    % get channel information
    [chan_info] = get_channel_info(strcat(sel_data.work_dir, ...
        '/mseed_files/', char(selected_stations(i)),'_',...
        char(stations.nc(st_ind)), '/channel_info.txt'),selected_stations(i));
    
    %read log file containing successfully downloaded data
    fileID2 = fopen([sel_data.work_dir,'/data_request_',...
        char(selected_stations(i)),'_',char(stations.nc(st_ind))]);
    C2 = textscan(...
        fileID2,'%f %f %f %f %f %f %f %f %f %f %f %f %s %s');
    fclose(fileID2);
    files = C2{1,14};
    
    % create waitbar
    wait_string = sprintf('pre-processing data for station %s ...',...
        char(selected_stations(i)));
    
    h = waitbar(0,wait_string,'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(h,'canceling',0);
    
    % loop around all files per station
    for n =  prev_n:length(files)
        
        if getappdata(h,'canceling')
            break
        end
        waitbar(n /length(files))
        
        % read event from text file & compare with channel info to find out
        % whether miniseed must be rotated because of misalignment
        [event] = read_event_file(C2,n);
        eventdatevec = datevec(...
            [num2str(event.year) '-' num2str(event.month) '-' num2str(event.day)]);
        eventnum = datenum(eventdatevec);
        if length(chan_info) > 1
            found_ch = 0;
            for j = 1:length(chan_info)
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
        file_mseed = extract_adv_sac(file_read);
        
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
                
                % rotate if nescessary
                if chan_info(choi).rot_flag
                    [event.n_amp,event.e_amp] = rot_az(event.n_amp,...
                        event.e_amp,chan_info(choi).cor_deg);
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
                        if  isfield(event.phases(i_phase),'cordeg')
                            to_save = to_save +1;
                            event.phases_to_analyze(to_save) = i_phase;
                        end
                    end
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
    delete(h);
    if exist('data','var') == 1
        save([dir_save,char(selected_stations(i)),'.mat'],'data');
        save([dir_save,char(selected_stations(i)),'_last_event.mat'],'n');
    else
        disp('there are no events for this station')
    end

    clear data n 
end

end