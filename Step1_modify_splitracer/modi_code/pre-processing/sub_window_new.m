function[ ] = sub_window_new(fig_id,an1,an2,an3,x,y1,y2,y3,...
    phases_t,sel_phase,hp)

% plot components for quality check

% usage:
% fig_id: number of plot
% an1,an2,an3: annotatoions for traces
% x: time vector
% y1,y2,y3: amplitudes
% phases_t: strcut with phases per event
% sel_phase: selected phase
% hp: parent figure container

% Copyright 2016 M.Reiss and G.Rümpker

% define size and position of plot window
seis_height = 0.8/7;
seis_width = 0.85/2;

seis_left = 0.05;
seis_bottom = 0.1;

% get subplot number which defines the trace number and the position of 
% the plot window 

no_seis = fig_id;
if fig_id ==1 || fig_id ==4 || fig_id ==7
    y=y1;
    an=an1;
elseif fig_id ==2 || fig_id ==5 || fig_id ==8
    y=y2;
    an=an2;
elseif fig_id ==3 || fig_id ==6 || fig_id ==9
    y=y3;
    an=an3;
end

% norm trace
N=length(x);

y1max=max(abs(y1));
y2max=max(abs(y2));
y3max=max(abs(y3));
ymax=max(y1max,y2max);
ymax=max(ymax,y3max);

xmax=max(abs(x));
xmin=min(abs(x));

y=y/ymax;

if fig_id < 4
    seis_bottom = 0.21;
end

% calculate position

if fig_id > 6 && fig_id < 10
    seis_left = seis_left*2+seis_width;
    seis_bottom = 0.21;
    no_seis= fig_id - 6;
end

%% plot trace

positionVector1 = [seis_left, seis_bottom + seis_height*(6-no_seis),...
    seis_width,seis_height];
axes('Parent',hp);
subplot('Position',positionVector1)

plot(x,y)
text(x(1)+0.05*(x(N)-x(1)),(1-0.5),an)

%% phase annotation
hold on

% only plot phases inside time window

time_min = min(x);
time_max = max(x);

for i_p = 1:length(phases_t)
    ttxks(i_p) = phases_t(i_p).tt_abs;
end
aa(1)=-1; aa(2)=1;
for it=1:length(ttxks)
    tt(1)=datenum(ttxks(it)); tt(2)=datenum(ttxks(it));
    if tt(1)>time_min && tt(1)<time_max
        plot(tt,aa,'-k')
        text(datenum(ttxks(it)),((-1)^it)*(1-0.14),char(phases_t(it).name))
    end
end

% plot selected phase in different color
tt1(1)=datenum(sel_phase(1)); tt1(2)=datenum(sel_phase(1));
plot(tt1,aa,'-g')
axis([xmin, xmax, -1, 1])

%% if time window was changed -  plot new window

if fig_id == 9 && length(sel_phase)>1
    
    aw1(1) = sel_phase(2);
    aw1(2) = sel_phase(2);
    aw2(1) = sel_phase(3);
    aw2(2) = sel_phase(3);
    plot(aw1,aa,'-r', aw2,aa,'-r')
    
elseif fig_id == 9
    sel_phase=datetime(datestr(sel_phase),'Format','d-MMM-y HH:mm:ss','Locale','en_US');
    aw1(1) = sel_phase - seconds(5);
    aw1(2) = sel_phase - seconds(5);
    aw2(1) = sel_phase + seconds(25);
    aw2(2) = sel_phase + seconds(25);
    plot(datenum(aw1),aa,'--r', datenum(aw2),aa,'--r')
end

hold off

%%  adjust axes

% adjust x axis
datetick('x',13)
axis([xmin, xmax, -1, 1])

% only plot x-axis annotation for every third plot
if fig_id ==1 || fig_id ==2 || fig_id ==4 || fig_id ==5 || fig_id ==7 ...
        || fig_id ==8
    set(gca,'XTickLabel','')
    axis([xmin, xmax, -1, 1])
end

end



