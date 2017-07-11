% Network 1DVAR+RTTOV retrieval: Plot netcdf output files
%
% Net1DPlot_netcdf plots data saved into netcdf output files 
%
% ncfile = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/NetCDF/level2/pay/tops_pay_mwrBL00_l2_hua_v00_20140101000353.nc';
% Net1DPlot_netcdf(ncfile)
%
% ncfile = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/NetCDF/level2/pay/tops_pay_mwrBL00_l2_ta_v00_20140101000353.nc';
% Net1DPlot_netcdf(ncfile)

function Net1DPlot_netcdf(ncfile)

MissingValue = -999; % should match C.MissingValue

% if contains(ncfile,'_l1_')
%    %plotlv1ncdata();
% end
% 
% if contains(ncfile,'_l2_ta_')
%    %plotlv2ncdata_ta();
% end

% Quick & dirty

% This plot ta
if contains(ncfile,'_l2_ta_')
    
indx = strfind(ncfile,'/');
ncfilename = ncfile(indx(end)+1:end);
station = ncfilename(6:8);

% load time
time = ncread(ncfile,'time'); % 'seconds since 1970-01-01 00:00:00 UTC'
[timestr,julday,datetime] = computertime(time);
hour = ( julday - floor(julday(1)) ) * 24;

% load profiles
z = ncread(ncfile,'height'); % altitude [m]
ta = ncread(ncfile,'ta');  % air_temperature [K]
ta_err = ncread(ncfile,'ta_err');  % air_temperature uncertainty [K]

% eliminate missing values
ta( find( ta == MissingValue ) ) = NaN;
ta_err( find( ta_err == MissingValue ) ) = NaN;

% plot (time series)
figure 
subplot(2,1,1)
pcolor(hour,z/1e3,ta); hc = colorbar;
xlabel('Time [h]'); ylabel('Height [km asl]'); set(get(hc,'Label'),'String','K','FontWeight','bold');
title(['1DVAR Ret Ta (' station ') ' datestr(datenum(datetime(1,1:3)))]);
shading flat; ylim([0 10]); xlim([0 24]); set(gca,'xtick',0:2:24);
subplot(2,1,2)
pcolor(hour,z/1e3,ta_err); hc = colorbar;
xlabel('Time [h]'); ylabel('Height [km asl]'); set(get(hc,'Label'),'String','K','FontWeight','bold');
title(['1DVAR Ret Ta uncertainty']);
shading flat; ylim([0 10]); xlim([0 24]); set(gca,'xtick',0:2:24);

% plot (profiles)
figure
ns = 4; % number of subplots
nr = length(ta(1,:));
np = floor( nr / ns );
for ip = 1:ns
    ir = ip * np;
    %subplot(2,ns/2,ip); plot(ta(:,ir),z/1e3,ta(:,ir)-ta_err(:,ir),z/1e3,':',ta(:,ir)+ta_err(:,ir),z/1e3,':');
    subplot(2,ns/2,ip); errorbar(ta(:,ir),z/1e3,ta_err(:,ir),'horizontal','CapSize',1); grid on;
    title(datestr(datenum(datetime(ir,1:6))));
    if ip==1 | ip==ns/2+1;
            ylabel('Height [km asl]');
    end
    if ip>=ns/2+1;
            xlabel('Ta [K]');
    end
end
suptitle(['1DVAR Ret Ta (' station ')']);

end

% This plot ha
if contains(ncfile,'_l2_hua_')
    
indx = strfind(ncfile,'/');
ncfilename = ncfile(indx(end)+1:end);
station = ncfilename(6:8);

% load time
time = ncread(ncfile,'time'); % 'seconds since 1970-01-01 00:00:00 UTC'
[timestr,julday,datetime] = computertime(time);
hour = ( julday - floor(julday(1)) ) * 24;

% load profiles
z = ncread(ncfile,'height'); % altitude [m]
hua = ncread(ncfile,'hua');  % air_temperature [K]
hua_err = ncread(ncfile,'hua_err');  % air_temperature uncertainty [K]

% eliminate missing values
hua( find( hua == MissingValue ) ) = NaN;
hua_err( find( hua_err == MissingValue ) ) = NaN;

% plot (time series)
figure 
subplot(2,1,1)
pcolor(hour,z/1e3,hua); hc = colorbar;
xlabel('Time [h]'); ylabel('Height [km asl]'); set(get(hc,'Label'),'String','kg/m^3','FontWeight','bold');
title(['1DVAR Ret HUa (' station ') ' datestr(datenum(datetime(1,1:3)))]);
shading flat; ylim([0 10]); xlim([0 24]); set(gca,'xtick',0:2:24);
subplot(2,1,2)
pcolor(hour,z/1e3,hua_err); hc = colorbar;
xlabel('Time [h]'); ylabel('Height [km asl]'); set(get(hc,'Label'),'String','kg/m^3','FontWeight','bold');
title(['1DVAR Ret Ta uncertainty']);
shading flat; ylim([0 10]); xlim([0 24]); set(gca,'xtick',0:2:24);

% plot (profiles)
figure
ns = 4; % number of subplots
nr = length(hua(1,:));
np = floor( nr / ns );
for ip = 1:ns
    ir = ip * np;
    %subplot(2,ns/2,ip); plot(ta(:,ir),z/1e3,ta(:,ir)-ta_err(:,ir),z/1e3,':',ta(:,ir)+ta_err(:,ir),z/1e3,':');
    subplot(2,ns/2,ip); errorbar(hua(:,ir),z/1e3,hua_err(:,ir),'horizontal','CapSize',1); grid on;
    title(datestr(datenum(datetime(ir,1:6))));
    if ip==1 | ip==ns/2+1;
            ylabel('Height [km asl]');
    end
    if ip>=ns/2+1;
            xlabel('HUa [kg/m^3]');
    end
end
suptitle(['1DVAR Ret HUa (' station ')']);
    
end


return