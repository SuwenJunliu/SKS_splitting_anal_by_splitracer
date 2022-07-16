function [trace]= extract_sac(file)

% this function reads the mutilple sac file, double entries
% a low pass of 1 s is administered to avoid aliasing effects
% that's only substitude the module for reading mseed


% usage:
% file: mseed file

% Copyright 2022 Junliu Suwen @ IGGCAS
postfix = '.sac';
comp_names = ["BHE","BHN","BHZ"];
sclice = 100;

try

    for i_comp = 1:length(comp_names)
        clear comp_name
        clear comp
        clear trace0
        clear dtfigure
        clear t
        clear DupIndex
        
        sac_file_name = [file '.' char(comp_names(i_comp)) postfix ];
        s = readsac(sac_file_name);
        
        if isfield(s,"KCMPNM")
            if s.KCMPNM ~= comp_names(i_comp)
                disp('the comp is not identical in sac head file?')
            end
            trace0 = s.DATA1;
           
            %demean
            trace0 = trace0-mean(trace0);
            % low pass butterworth filter
            dt = s.DELTA;
            
            
            
            %trace0 = buttern_low(trace0',6,1,dt);
            %trace0 = buttern_low(trace0',6,1,dt);
            trace0 = buttern_filter(trace0',6,0.01,1,dt);
            
            
            [year,month,day] = cal_day(s.NZYEAR,s.NZJDAY);
            starttime = datenum(year,month,day,s.NZHOUR,s.NZMIN,s.NZSEC + s.B);
            endtime =  datenum(year,month,day,s.NZHOUR,s.NZMIN,s.NZSEC + s.E); 
            
            
            t = linspace(starttime,endtime,s.NPTS);
        end
        
        if strcmp(comp_names(i_comp),'BHE')
            trace(1).comp = 'East';
            trace(1).amp = trace0;
            trace(1).time = t;
            trace(1).dt = dt;
        end
        
        if strcmp(comp_names(i_comp),'BHN')
            trace(2).comp = 'North';
            trace(2).amp = trace0;
            trace(2).time = t;
            trace(2).dt = dt;
        end
        if strcmp(comp_names(i_comp),'BHZ')
            trace(3).comp = 'Z';
            trace(3).amp = trace0;
            trace(3).time = t;
            trace(3).dt = dt;
        end
        
    end
    
    return

catch
    disp('could not read sac, check the file name or "postfix" in func')
    trace = 0;
    return
end





end