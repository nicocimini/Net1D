% O minus B monitoring: plot stats
%
% OmBPlotStats plots the O-B stats as compute by OmBComputeStats.

function OmBPlotStats(OmB,OmBpath)

savefig = 1; grey = [0.5 0.5 0.5];
clrs = ['k', 'b', 'c', 'r', 'g', 'm'];
na = length(OmB.angles_el);
nc = length(OmB.channels);
switch OmB.instrument 
    case 'HATPRO'
       kband = 1:7;
       vband = 8:14;
    case 'MP3000A'
       kband = 1:5;
       vband = 6:12;
end
    
% set out folder
if savefig
   fig_path = [OmBpath 'FIGS/']; 
   if ~exist(fig_path,'dir'); mkdir(fig_path); end;
end


% time 
[timestr,julday] = computertime(OmB.time);
timestr = datestr([timestr(1,:); timestr(end,:)],'YYYYmmddhh');

% cloudy and clear sky
cloudy = find(OmB.cld31);
clrsky = find(OmB.cld31==0);
zen_only = 1;

% Plot Stats
% All channels
figure
for ia = 1:na
    subplot(na,1,ia)
    plot(OmB.channels,OmB.mean(:,ia),'k*-',OmB.channels,OmB.stdv(:,ia),'ro-'); grid on;
    if ia==1; title(['O-B stats for ' OmB.station_id '(' OmB.instrument ') from ' timestr(1,:) ' to ' timestr(2,:)]); end;
    if ia==na; xlabel('Channels [GHz]'); end;
    if ia==na/2; ylabel('Mean and STD Tb diff [K]'); end
    text(60.5,0,[num2str(OmB.angles_el(ia),'%5.2f') '\circ']);
end
if savefig; saveas(gcf,[fig_path 'OmB_Stats_' OmB.station_id '_from_' timestr(1,:)  '_to_' timestr(2,:)],'fig'); end;

% K and V channels separated
figure
for ia = 1:na
    subplot(na,2,2*ia-1)
    plot(OmB.channels(kband),OmB.mean(kband,ia),'k*-',OmB.channels(kband),OmB.stdv(kband,ia),'ro-'); grid on;
    if ia==1; title(['K-band' ]); end;
    if ia==na; xlabel('Channels [GHz]'); end;
    if ia==na/2; ylabel('Mean and STD Tb diff [K]'); end
    text(60.5,0,[num2str(OmB.angles_el(ia),'%5.2f') '\circ']);
    subplot(na,2,2*ia)
    plot(OmB.channels(vband),OmB.mean(vband,ia),'k*-',OmB.channels(vband),OmB.stdv(vband,ia),'ro-'); grid on;
    if ia==1; title(['V-band']); end;
    if ia==na; xlabel('Channels [GHz]'); end;
    text(60.5,0,[num2str(OmB.angles_el(ia),'%5.2f') '\circ']);
end
suptitle(['O-B stats for ' OmB.station_id '(' OmB.instrument ') from ' timestr(1,:) ' to ' timestr(2,:)]);
if savefig; saveas(gcf,[fig_path 'OmB_Stats_' OmB.station_id '_from_' timestr(1,:)  '_to_' timestr(2,:) 'KandVchan'],'fig'); end;

% Plot Time Series
figure
diff = OmB.diff(:,zen_only,:);
for ic = 1:nc
    subplot(nc,1,ic)
    %plot(julday,diff(ic,:),'.','Color',grey); grid on; hold on;
    plot(julday(clrsky),diff(ic,clrsky),'k.'); grid on; hold on;
    if ic==1; title(['O-B time series for ' OmB.station_id '(' OmB.instrument ') from ' timestr(1,:) ' to ' timestr(2,:)]); end;
    if ic==nc; xlabel('Julian day [d]'); end;
    if ic==nc/2; ylabel('Tb diff [K]'); end
    text(julday(end)+0.2,0,[num2str(OmB.channels(ic),'%5.2f') 'GHz']);
end
if savefig; saveas(gcf,[fig_path 'OmB_TimeS_' OmB.station_id '_from_' timestr(1,:)  '_to_' timestr(2,:)],'fig'); end;

% Plot Stats as Pauline
figure
subplot(2,1,1)
   plot(OmB.channels(kband),OmB.mean(kband,1),'kx-',OmB.channels(kband),OmB.stdv(kband,1),'ro--'); grid on;
   title(['O-B stats for ' OmB.station_id '(' OmB.instrument ') from ' timestr(1,:) ' to ' timestr(2,:)]);
   legend(['Mean @ ' num2str(OmB.angles_el(1),'%5.2f')],['STD @ ' num2str(OmB.angles_el(1),'%5.2f')]);
   xlabel('Channels [GHz]'); ylabel('Mean and STD Tb O-B diff [K]'); xlim([22 32])
subplot(2,1,2)
   for ia = 1:na % mean
       plot(OmB.channels(vband),OmB.mean(vband,ia),[clrs(ia) 'x-']); hold on;
   end
   for ia = 1:na % stdv
       plot(OmB.channels(vband),OmB.stdv(vband,ia),[clrs(ia) 'o--']); hold on;
   end
   legend(num2str(OmB.angles_el(:),'%5.2f'),'Location','southeast'); grid on;
   xlabel('Channels [GHz]'); ylabel('Mean and STD Tb O-B diff [K]'); xlim([51 59]);
if savefig; saveas(gcf,[fig_path 'OmB_Stats_asPauline' OmB.station_id '_from_' timestr(1,:)  '_to_' timestr(2,:)],'png'); end;

% Plot Time Series as Pauline
figure
chan_sel = [1 kband(end)-1 kband(end) vband(2) vband(end)];
nc = length(chan_sel);
for ic = 1:nc
    subplot(nc,1,ic)
    plot(julday(clrsky),diff(chan_sel(ic),clrsky),'k*-'); grid on; hold on; ylim([-10 10]); set(gca,'ytick',[-10:5:10]);
    %legend(num2str(OmB.channels(chan_sel(ic)),'%5.2f')); ylabel('O-B Tb [K]'); ylim([-10 10]);
    text(julday(end)+0.2,0,[num2str(OmB.channels(chan_sel(ic)),'%5.2f') 'GHz']); ylabel('O-B Tb [K]'); 
    if ic==1; title(['O-B time series for ' OmB.station_id '(' OmB.instrument ') from ' timestr(1,:) ' to ' timestr(2,:)]); end;
    if ic==nc; xlabel('Julian day [d]'); end;
end
if savefig; saveas(gcf,[fig_path 'OmB_TimeS_asPauline' OmB.station_id '_from_' timestr(1,:)  '_to_' timestr(2,:)],'png'); end;

return