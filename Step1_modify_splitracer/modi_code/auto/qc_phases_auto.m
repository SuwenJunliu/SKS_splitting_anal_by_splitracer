function [saved_events] = qc_phases_auto(data,prev_events,...
    start_index,sel_data)

% quality check for phases after inital preparation and automated window
% selection
% usage:
% data: saved events/phases from previous preparation step
% station: station name
% snr_cutoff: previously defined snr cut off

% Copyright 2021 F.Link, M.Reiss and G.Rümpker 

% read variable names
events = fieldnames(data);

to_keep = 0;

axidopts = 'abcdefghijklmnopqrstuvwxyz';

% check with which event to start if data was added later & previous
% analysis exists 

if isstruct(prev_events)
    final_events = prev_events;
     to_keep = sum(cell2mat(strfind(fieldnames(prev_events),'event')));
end

% loop over all events in data
singleflag = 0; % change to 1 if only 1 phase per event should be used
for iEvents = start_index:length(events)

    % choose event to be analyzed
    fn = events{iEvents};
    event = data.(fn);
    phases_t = event.phases;

    % choose phases to be analyzed
    Nph = length(event.phases_to_analyze);
    for an_phases = 1:Nph
        
        naxiopt = 1;

        i_phase = event.phases_to_analyze(an_phases);

        % cut timerange
        t_xks = event.phases(i_phase).tt_abs;
        

        cutc = datenum((t_xks-seconds(50)));
        [~, index] = min(abs(event.time-cutc));
        

        index_end = index+(100/data.(fn).sr);
        cut.time = event.time(index:index_end);
        
        bhw = tukeywin(length(index:index_end),0.3);
        bhw = ones(size(bhw));
        cut.z_unfilt = event.z_amp(index:index_end).*bhw';
        cut.north_unfilt = event.n_amp(index:index_end).*bhw';
        cut.east_unfilt = event.e_amp(index:index_end).*bhw';
        
        % check if a channel is dead
        if (max(cut.north_unfilt) == 0 && min(cut.north_unfilt) == 0) || ...
                (max(cut.east_unfilt) == 0 && min(cut.east_unfilt) == 0) || ...
                (max(cut.z_unfilt) == 0 && min(cut.z_unfilt) == 0)
            disp('Traces with dead channel. Skipping event.')
            continue
        end

        % filter according to filter range
        zfilttemp = buttern_filter(event.z_amp,2,...
            1/sel_data.p2,1/sel_data.p1,event.sr);
        nfilttemp = buttern_filter(event.n_amp,2,...
            1/sel_data.p2,1/sel_data.p1,event.sr);
        efilttemp = buttern_filter(event.e_amp,2,...
            1/sel_data.p2,1/sel_data.p1,event.sr);
        cut.z = zfilttemp(index:index_end).*bhw;
        cut.north = nfilttemp(index:index_end).*bhw;
        cut.east = efilttemp(index:index_end).*bhw;
        
        % filter in long period range
        cut.north_lp = buttern_filter(cut.north_unfilt,2,...
            1/50,1/15,event.sr);
        cut.east_lp = buttern_filter(cut.east_unfilt,2,...
            1/50,1/15,event.sr);

        % rotate to radial/transverse
        [cut.radc, cut.trac] = rad_tra(cut.north, cut.east,event.baz);
        fnorm = max(cut.radc);
        cut.radc = cut.radc/fnorm;
        cut.trac = cut.trac/fnorm;

        % read phases and prepare maximum time range
        x = cut.time;

        time_min = min(x);
        time_max = max(x);
        maxtimeph = time_max;
        
        clear ttxks phindx
        for i_p = 1:length(phases_t)
            ttxks(i_p) = phases_t(i_p).tt_abs;
        end
        for it=1:length(ttxks)
            tt(1)=datenum(ttxks(it)); tt(2)=datenum(ttxks(it));
            if tt(1)>time_min && tt(1)<time_max
                if strcmp(char(phases_t(it).name),'S')
                    maxtimeph = min([maxtimeph tt(1)]);
                end
                if strcmp(char(phases_t(it).name),'ScS')
                    maxtimeph = min([maxtimeph tt(1)]);
                end
                if strcmp(char(phases_t(it).name),'Sdiff')
                    maxtimeph = min([maxtimeph tt(1)]);
                end
            end
            if min(cut.time) < datenum(ttxks(it)) && max(cut.time) > datenum(ttxks(it)) 
            phindx(it) = floor(interp1(cut.time,1:length(cut.time),datenum(ttxks(it))));
            else
                phindx(it) = 0;
            end
        end
        maxindxph = floor(interp1(cut.time,1:length(cut.time),maxtimeph));

        % prepare time frequency analysis
        T = event.sr;
        Fs = 1./T;
        Fcutoff = 1./sel_data.p1;
        Fcutofflong = 1./sel_data.p2;
        
        endmaxindx = min([maxindxph length(cut.radc)]);
        
        Fsteps = 100;

        ff = linspace(Fcutofflong,Fcutoff*2,Fsteps);
        maxwin = round(1/max(ff)*Fs*1.5);
        S = zeros(Fsteps-1,endmaxindx-(maxwin-1));
        FF = zeros(1,Fsteps);
        for i = 1:(length(ff)-1)
            ffoi = ff(i:i+1);
            Ntkw = round(1/max(ffoi)*Fs*1.5);
            bhw = tukeywin(Ntkw,0.3);
            if Ntkw > endmaxindx
                S(i,:) = 0;
                TTtemp = 0;
                F = ffoi(1);
            else
                [ss,F,TTtemp]=spectrogram(cut.radc(1:endmaxindx),bhw,Ntkw-1,ffoi,Fs);
                [ssT,~,~]=spectrogram(cut.trac(1:endmaxindx),bhw,Ntkw-1,ffoi,Fs);
                S(i,(round(Ntkw/2)-round(maxwin/2)+(1:(length(ss)))))= ss(1,:)/Ntkw+ssT(1,:)/Ntkw;
            end
            FF(i) = F(1);
            if i == 1
                TT = TTtemp;
            else
                if length(TT) < length(TTtemp)
                    TT = TTtemp;
                end
            end
        end
        S(length(ff),(round(Ntkw/2)-round(maxwin/2)+(1:(length(ss)))))= ss(2,:)/Ntkw;
        FF(length(ff)) = F(2);
        [~,position]=max(abs(S(:)));
        [I,J]=ind2sub(size(S),position);
        SS = abs(S);
        
        %qf = 0.8;
        %qt = 0.5;
        %modi 
        qf = 0.8;
        qt = 0.65;
        
        Smax = max(SS,[],2);
        indxoi = 1:length(Smax);
        Soi = sum(SS(Smax>max(SS(:))*qf,:),1);
        Soi = Soi./max(Soi(:));
        Fb = FF(Smax>max(SS(:))*qf);
        Ioi = indxoi(Smax>max(SS(:))*qf);
        II = Ioi(Fb==min(Fb));
        minI = min(indxoi(Smax>max(SS(:))*qf));
        maxI = max(indxoi(Smax>max(SS(:))*qf));
        [maxSoi,JJ] = max(Soi);
        indx = 1:length(Soi);
        indxoi = indx(Soi>=maxSoi*qt);
        indx2 = 1:length(indxoi);
        gaps = indx2(abs(indxoi(1:end-1)-indxoi(2:end))>1);
        gapss = indxoi(gaps+1);
        gapsb = indxoi(gaps);
        gapsoi = gaps(gapss<JJ);
        gapss = gapss(gapss<JJ);
        gapss((gapss-indxoi(gapsoi))<200) = [];
        gapsoi = gaps(gapsb>JJ);
        gapsb = gapsb(gapsb>JJ);
        gapsb((indxoi(gapsoi+1)-gapsb)<50) = [];
        if ~isempty(gapss)
            resstartindx = max(gapss)-round((5+1/FF(II))*Fs/2)+round(maxwin/2);
            resstartindx1 = max(gapss)+round(maxwin/2);
        else
            resstartindx = min(indxoi)-round((5+1/FF(II))*Fs/2)+round(maxwin/2);
            resstartindx1 = min(indxoi)+round(maxwin/2);
        end
        if ~isempty(gapsb)
            resendindx = min(gapsb)+round((5+1/FF(II))*Fs/2)-1+round(maxwin/2);
            resendindx1 = min(gapsb)+round(maxwin/2);
        else
            resendindx = max(indxoi)+round((5+1/FF(II))*Fs/2)-1+round(maxwin/2);
            resendindx1 = max(indxoi)+round(maxwin/2);
        end
        resstartindx = max([1 resstartindx]);
        resendindx = min([resendindx length(cut.radc)]);
        
        FoSi = 8;
        FoSiTi = 10;
        fig = figure('Visible','off','Position',[50 50 600 700]);
        sp1 = subplot(5,1,1);
        time = 0:event.sr:event.sr*length(cut.radc);
        plot(time(time>min(TT)&time<max(TT)),cut.radc(time>min(TT)&time<max(TT)),'b');
        hold on
        plot([time(resstartindx) time(resstartindx)],[min(cut.radc) max(cut.radc)],'r')
        plot([time(min([length(time) resendindx])) time(min([length(time) resendindx]))],[min(cut.radc) max(cut.radc)],'r')
        for i = 1:length(phindx)
            if phindx(i) ~= 0
                plot([time(phindx(i)) time(phindx(i))],[min(cut.radc) max(cut.radc)],'g');
                text(time(phindx(i))+1,0.5*max(cut.radc),char(phases_t(i).name),'FontSize',FoSi)
            end
        end
        axis([min(TT) max(TT) min(cut.radc) max(cut.radc)]) 
        ylabel('Amplitude','FontSize',FoSi)
        title('Radial Component','FontSize',FoSiTi)
        gp = get(gca,'position');
        dim = [gp(1)-0.05 gp(2)+gp(4)*1.3 0 0];
        str = [axidopts(naxiopt) ')'];
        naxiopt = naxiopt+1;
        annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none','FontSize',FoSiTi,'FontWeight','bold');
        
        sp2 = subplot(5,1,2);
        time = 0:event.sr:event.sr*length(cut.radc);
        plot(time(time>min(TT)&time<max(TT)),cut.trac(time>min(TT)&time<max(TT)),'b');
        hold on
        plot([time(resstartindx) time(resstartindx)],[min(cut.radc) max(cut.radc)],'r')
        plot([time(min([length(time) resendindx])) time(min([length(time) resendindx]))],[min(cut.radc) max(cut.radc)],'r')
        axis([min(TT) max(TT) min(cut.radc) max(cut.radc)]) 
        ylabel('Amplitude','FontSize',FoSi)
        title('Transverse Component','FontSize',FoSiTi)
        gp = get(gca,'position');
        dim = [gp(1)-0.05 gp(2)+gp(4)*1.3 0 0];
        str = [axidopts(naxiopt) ')'];
        naxiopt = naxiopt+1;
        annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none','FontSize',FoSiTi,'FontWeight','bold');
        
        sp3 = subplot(5,1,3); 
        contourf(TT,log10(FF),abs(S),20,'LineStyle','none')
        hold on
        plot([TT(1) TT(end)],[log10(FF(I)) log10(FF(I))],'r')
        plot([TT(1) TT(end)],[log10(FF(minI)) log10(FF(minI))],'g-.')
        plot([TT(1) TT(end)],[log10(FF(maxI)) log10(FF(maxI))],'g-.')
        plot([time(resstartindx1) time(resstartindx1)],[min(log10(FF)) max(log10(FF))],'g')
        plot([time(resendindx1) time(resendindx1)],[min(log10(FF)) max(log10(FF))],'g')
        plot([time(resstartindx) time(resstartindx)],[min(log10(FF)) max(log10(FF))],'r')
        plot([time(min([length(time) resendindx])) time(min([length(time) resendindx]))],[min(log10(FF)) max(log10(FF))],'r')
        axis([min(TT) max(TT) log10(FF(1)) log10(FF(end))])
        ylabel('Frequency log(F)]','FontSize',FoSi)
        title('Spectral Energy','FontSize',FoSiTi)
        gp = get(gca,'position');
        dim = [gp(1)-0.05 gp(2)+gp(4)*1.3 0 0];
        str = [axidopts(naxiopt) ')'];
        naxiopt = naxiopt+1;
        annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none','FontSize',FoSiTi,'FontWeight','bold');
        
        sp4 = subplot(5,1,4);
        plot(TT,Soi,'b')
        hold on
        plot([time(1) time(end)],[qt*max(Soi(:)) qt*max(Soi(:))],'r')
        plot([time(resstartindx1) time(resstartindx1)],[min(Soi) max(Soi)],'g')
        plot([time(resendindx1) time(resendindx1)],[min(Soi) max(Soi)],'g')
        plot([time(resstartindx) time(resstartindx)],[min(Soi) max(Soi)],'r')
        plot([time(min([length(time) resendindx])) time(min([length(time) resendindx]))],[min(Soi) max(Soi)],'r')
        axis([min(TT) max(TT) min(Soi) max(Soi)])
        xlabel('Time in seconds [s]','FontSize',FoSi)
        ylabel('Spectral Energy','FontSize',FoSi)
        title('Spectral Energy for reduced Frequencies','FontSize',FoSiTi)
        gp = get(gca,'position');
        dim = [gp(1)-0.05 gp(2)+gp(4)*1.3 0 0];
        str = [axidopts(naxiopt) ')'];
        naxiopt = naxiopt+1;
        annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none','FontSize',FoSiTi,'FontWeight','bold');
        
        filepath = [sel_data.work_dir '/graphics_output/pre-processing/Filter_' num2str(sel_data.p1) '-' num2str(sel_data.p2) 's_SNR_' num2str(sel_data.snr) '/' char(sel_data.station{1}) '/Spectrograms/'];
        if ~exist(filepath,'dir')
            mkdir(filepath)
        end
        filename = [filepath 'event' num2str(iEvents) '_phase' char(phases_t(i_phase).name)];
        
        
        % Check noise level
        SNRflag = 1;
        trace1 = cut.radc(max([1 resstartindx]):min([end resendindx]));
        trace2b = cut.radc(1:min([(resstartindx-1) (length(trace1)-1)]));
        maxlev = sum(trace1.^2)/length(trace1);
        minlev = sum(trace2b.^2)/length(trace2b);

        % Check position
        keepitflag = 1;
        edgeflag = 0;
        text(min(TT)+1,0.8,datestr(event.origin_time,'yyyy-mm-dd HH:MM:SS'),'Color','Red','FontSize',FoSi)
        if resendindx >= maxindxph || resendindx >= endmaxindx
            disp('XKS-phase at edge of realistic time range, Skipping trace. (End of trace)')
            edgeflag = 1;
            subplot(sp1);
            text(min(TT)+1,0.6,'Time Range Error','Color','Red','FontSize',FoSi)
            keepitflag = 0;
        end
        if J == 1 && keepitflag%&& newI == 1
            disp('XKS-phase at edge of realistic time range, Skipping trace. (Start of trace)')
            edgeflag = 1;
            subplot(sp1);
            text(min(TT)+1,0.6,'Time Range Error','Color','Red','FontSize',FoSi)
            keepitflag = 0;
        end 
        if FF(I) < 1.5*Fcutofflong && keepitflag
            disp('XKS-phase at edge of realistic frequency range, Skipping trace')
            edgeflag = 1;
            subplot(sp1);
            text(min(TT)+1,0.6,'Frequency Range Error','Color','Red','FontSize',FoSi)
            keepitflag = 0;
        end
        if SS(II,resendindx1-round(maxwin/2)) == 0 && keepitflag
            disp('XKS-phase at edge of realistic time range, Skipping trace. (Maximum at end of trace)')
            edgeflag = 1;
            subplot(sp1);
            text(min(TT)+1,0.6,'Time Range Error','Color','Red','FontSize',FoSi);
            keepitflag = 0;
        end
        clear indxoi stalta sta lta
        for i = 1:(endmaxindx-round(20/event.sr))
            indxoi(i) = round(mean(i:(round(20/event.sr)+i-1)));
            sta(i) = mean(cut.radc(i:(round(20/event.sr)+i-1)).^2);
            if i < round(25/event.sr)
                lta(i) = mean(cut.radc(1:round(50/event.sr)).^2);
            elseif i > endmaxindx-round(25/event.sr)
                lta(i) = mean(cut.radc(endmaxindx-round(50/event.sr):endmaxindx).^2);
            else
                lta(i) = mean(cut.radc(i-round(25/event.sr)+1:i+round(25/event.sr)).^2);
            end
            stalta(i) = sta(i)/lta(i); 
        end
        [slmax,slItemp]=max(stalta(indxoi<resendindx&indxoi>resstartindx)); 
        if isempty(slmax)
            [slmax,slI] = max(stalta);
        else
            slI = max(1,slItemp+resstartindx-indxoi(1)+1);
        end
        
        sp5 = subplot(5,1,5);
        plot(TT(indxoi),stalta,'b')
        hold on;
        plot([TT(1) TT(indxoi(1))],[stalta(1) stalta(1)],'b--')
        plot([TT(indxoi(end)) TT(end)],[stalta(end) stalta(end)],'b--')
        plot([TT(indxoi(slI)) TT(indxoi(slI))],[min(stalta)*0.9 max(stalta)*1.1],'g')
        plot([min(TT) max(TT)],[2.1 2.1],'r')
        axis([min(TT) max(TT) min(stalta)*0.9 max(stalta)*1.1])
        xlabel('Time in seconds [s]','FontSize',FoSi)
        ylabel('STA/LTA','FontSize',FoSi)
        title('Ratio of STA over LTA','FontSize',FoSiTi)
        gp = get(gca,'position');
        dim = [gp(1)-0.05 gp(2)+gp(4)*1.3 0 0];
        str = [axidopts(naxiopt) ')'];
        naxiopt = naxiopt+1;
        annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none','FontSize',FoSiTi,'FontWeight','bold');
        staltaflag = 1;
        if slmax < 2 && keepitflag
            %modi slmax from 2.1 to 2
            disp('XKS-phase shows low STA/LTA, Skipping trace.')
            staltaflag = 0;
            subplot(sp1);
            text(min(TT)+1,0.6,'STA/LTA Error','Color','Red','FontSize',FoSi);
            keepitflag = 0;
        end
        if (indxoi(slI) > resendindx || indxoi(slI) < resstartindx) && keepitflag
            disp('XKS-phase shows low STA/LTA, Skipping trace.')
            staltaflag = 0;
            subplot(sp1);
            text(min(TT)+1,0.6,'STA/LTA Error','Color','Red','FontSize',FoSi);
            keepitflag = 0;
        end
        if ~edgeflag && SNRflag && staltaflag
            disp('Good XKS-phase')
            subplot(sp1);
            text(min(TT)+1,0.6,'Accepted','Color','Green','FontSize',FoSi)
            event.phases(i_phase).tw(1) = datenum(cut.time(resstartindx));
            event.phases(i_phase).tw(2) = datenum(cut.time(resendindx));

            win = tukeywin(resendindx-resstartindx+1,0.3);
            cut.pm_east_sp = cut.east(resstartindx:resendindx).*win;
            cut.pm_north_sp = cut.north(resstartindx:resendindx).*win;
            cut.pm_east_lp = cut.east_lp(resstartindx:resendindx).*win;
            cut.pm_north_lp = cut.north_lp(resstartindx:resendindx).*win;
            ntemp_sp = zeros(size(cut.north));
            etemp_sp = ntemp_sp;
            ntemp_sp(1:length(cut.pm_north_sp)) = cut.pm_north_sp;
            etemp_sp(1:length(cut.pm_east_sp)) = cut.pm_east_sp;
            ntemp_pm = buttern_filter(cut.pm_north_sp,2,1/((resendindx-resstartindx)*event.sr),FF(minI),event.sr);
            etemp_pm = buttern_filter(cut.pm_east_sp,2,1/((resendindx-resstartindx)*event.sr),FF(minI),event.sr);
            fnorm = max(sqrt((ntemp_pm.^2)+(etemp_pm.^2)));
            ntemp_pm = ntemp_pm/fnorm;
            etemp_pm = etemp_pm/fnorm;
            cut.pm_north_lp = ntemp_pm;
            cut.pm_east_lp = etemp_pm;

            % calculate new short/long axis ratio
            [~,xlam1,xlam2] = covar(cut.pm_north_sp,cut.pm_east_sp);
            event.phases(i_phase).xlam_ratio = xlam2/xlam1;

            [baz_long]=pm_auto(cut.pm_north_lp,cut.pm_east_lp);

            % calculate new difference to backazimuth
            if event.baz > baz_long
                event.phases(i_phase).cordeg = mod(event.baz - baz_long,180);
            elseif event.baz <= baz_long
                event.phases(i_phase).cordeg = mod(event.baz - baz_long,-180);
            end

            if event.phases(i_phase).cordeg > 90
                event.phases(i_phase).cordeg = ...
                    event.phases(i_phase).cordeg-180;
            end

            if event.phases(i_phase).cordeg < -90
                event.phases(i_phase).cordeg = ...
                    event.phases(i_phase).cordeg +180;
            end

            % calculate new short/long axis ratio
            event.phases(i_phase).xlam_ratio = xlam2/xlam1;

            if to_keep == 0
                to_keep = to_keep +1;
                new_name = ['event' int2str(to_keep)];
                final_events.(new_name)=event;
            else
                fn = ['event' int2str(to_keep)];
                if final_events.(fn).origin_time == event.origin_time
                    final_events.(fn) = event;
                else
                    to_keep = to_keep +1;
                    new_name = ['event' int2str(to_keep)];
                    final_events.(new_name) = event;
                end
            end
        end
        
        print(fig,filename,'-r300','-dpng')
        close(fig)

        if singleflag
            return
        end
    end
end
if exist('final_events','var')
saved_events = final_events;
else
    saved_events = 0;
end
end