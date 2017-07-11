% This function may go in Net1DLoad_level1.m once it settles
% Or move Net1DLoad_level1_BL out of Net1DLoad_level1.m
% 
% lv1file = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/level1/ces/tops_ces_mwr00_l1_tb_v00_20140101000001.nc'
% sups_joy_mwr00_l1_tb_v01_20140101000004.nc
% sups_lac_mwr00_l1_tb_v00_20140101000218.nc
% tops_pay_mwr00_l1_tb_v00_20140101000001.nc
% sups_sir_mwr00_l1_tb_v01_20140102000000.nc
% tops_rao_mwr00_l1_tb_v01_20140101000000.nc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HATPRO and MP3000A Zenith-only files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [O,C] = Net1DLoad_level1_ZH(C,lv1file);

% update C
C.lat = ncread(lv1file,'lat'); % 'degree_north'
C.lon = ncread(lv1file,'lon'); % 'degree_east'
C.asl = ncread(lv1file,'zsl'); % 'altitude above mean sea level [m]'
% Global Attributes:
O.GA.Title = ncreadatt(lv1file,'/','Title');
O.GA.Institution = ncreadatt(lv1file,'/','Institution');
O.GA.Contact_person = ncreadatt(lv1file,'/','Contact_person');
O.GA.Source = ncreadatt(lv1file,'/','Source');
O.GA.History = ncreadatt(lv1file,'/','History');
O.GA.Conventions = ncreadatt(lv1file,'/','Conventions');
O.GA.Processing_date = ncreadatt(lv1file,'/','Processing_date');
O.GA.Author = ncreadatt(lv1file,'/','Author');
try O.GA.Comments = ncreadatt(lv1file,'/','Comments'); catch; O.GA.Comments = ' '; end;
try O.GA.License = ncreadatt(lv1file,'/','License'); catch; O.GA.License = ncreadatt(lv1file,'/','Licence'); end;
try O.GA.Measurement_site = ncreadatt(lv1file,'/','Measurement_site'); catch; O.GA.Measurement_site = upper(C.station); end;
% NC info:
O.ncinfo = ncinfo(lv1file);


% Find flagged data
% for Joyce: 'Flags indicate data that the user should only use with care. In cases of doubt, please refer to the contact person. A Fillvalue of 0 means that data has not been flagged. Bands refer to the measurement ranges (if applicable) of the microwave radiometer; i.e band 1: 20-30 GHz, band 2: 50-60 GHz, band 3: 90 GHz; tb valid range: [  2.70, 330.00] in K; '
% for Payerne: 'rain_flag radiometer_alarm; Discard data with rain_flag=1 and radiometer_alarm=1'
flag = ncread(lv1file,'flag');  
switch C.station_id
    case 'ces'
        badindx = find(flag);
    case 'joy'
        badindx = find(not(isnan(flag)));
    case 'lac'
        badindx = find(not(isnan(flag)));
    case 'pay'
        badindx = find(flag == 1);
    case 'rao'
        % For zh files, it's OK to get badindx from original flag (since ele is always 90)
        badindx = find(flag);
    case 'sir'
        badindx = find(not(isnan(flag)));
end

% load data
time = ncread(lv1file,'time'); % 'seconds since 1970-01-01 00:00:00 UTC'
azi = ncread(lv1file,'azi'); % '0=North, 90=East, 180=South, 270=West'
ele = ncread(lv1file,'ele'); % 'retrieval elevation angle'
freq = ncread(lv1file,'freq_sb'); % 'frequency of microwave channels [GHz]'
tbs = ncread(lv1file,'tb');  % 'brightness temperature [K]' - n_freq,time
ta = ncread(lv1file,'ta');  % 'air_temperature [K]'
pa = ncread(lv1file,'pa');  % 'air_pressure [Pa]'
hur = ncread(lv1file,'hur');  % 'relative_humidity [%/100]'
% comment - Some types of MWR require a systemmatic adjustement of the measured TBs. This variable gives the offset which was subtracted from each measurement. The offset was deteremined from COSMO-DE analysis.
offset_tb = ncread(lv1file,'offset_tb'); % 'brightness temperature offset subtracted from measured brightness temperature' [K]
% comment - This variable is an estimate of the one-standard-deviation calibration error to be expected from an absolute system calibration, i.e. the likely systematic error of brightness temperature. As a reference see Maschwitz et al. 2012, AMT (Tab. 5). However, these numbers differ from instrument to instrument and should be adapted accordingly. Values only valid for elevation angles larger than 20deg.
tb_bias = ncread(lv1file,'tb_bias'); % 'brightness temperature bias subtracted from measured brightness temperature' [K]
% comment - RPG offers a frequency shift within the radiometer software. This shift will modify the TBs calculated. Original TBs cannot be reconstructed. The variable given here is intended to inform the user which frequency shifts were applied to the given TBs
freq_shift = ncread(lv1file,'freq_shift'); % 'frequency shift applied to correct measured brightness temperature [GHz]'
% comment - This variable is calculated from brightness temperature observations of an internal black body whose physical temperature is known. The square root of the matrix diagonal gives the brightness temperature random error of each frequency channel. Values only valid for elevation angles larger than 20deg.
tb_cov = ncread(lv1file,'tb_cov'); % 'error covariance matrix of brightness temperature channels' [K^2]

% Should be zenith only
elangs = unique(ele);
azi = unique(azi);
if length(elangs) ~= 1; disp('something is wrong with the zenith pointing'); return; end;
ele = elangs;

% This only applies to MP3000A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(C.instrument,'MP3000A')
    % load extra fields
    wl_irp = ncread(lv1file,'wl_irp'); % 'infrared pyrometer central wavelength [micrometers]'
    tb_irp = ncread(lv1file,'tb_irp'); % 'infrared pyrometer brightness temperature [K]'
    el_irp = ncread(lv1file,'ele_irp'); % 'infrared pyrometer elevation angle [degree]'
    tb_irp = tb_irp'; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Maybe the following is not needed for zenith obs!
% Other instrument specific fix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To be removed one day we make the provider fix these %%%%%%%%%%%%%%%%%%%%
switch C.station_id
    case 'ces'
          %tbs = permute(tbs,[2 1 3]); % original files switch the nchn and nang dimensions
    case 'joy'

    case 'lac'

    case 'pay'

    case 'rao'

    case 'sir'
          %tbs = tbs(:,1:6,:); % discarding the observations at 0 elevation 
          %ele(1:6) = ele(1:6) - 9*1e4; % don't know the reason...
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discard flagged data
time(badindx) = []; 
tbs(:,badindx) = []; % just n_freq,time
ta(badindx) = [];
pa(badindx) = [];
hur(badindx) = [];
offset_tb(:,badindx) = []; % just n_freq,time
if strcmp(C.instrument,'MP3000A')
    tb_irp(badindx) = [];
    el_irp(badindx) = [];
end

% Pauline may want to have the original resolution to compute a cloud flag
% from Tb31, so the "Select time" section may be not wanted 
% Alternatively, we could leave as is, but compute here a cloud flag at reduced resolution
% based on Tb31 at original resolution (similarly within Net1DLoad_level1_BL)
% Discuss with Pauline

% Select time
% MAYBE THIS IS NOT NEEDED for ZH
[timestr,julday,datetime] = computertime(time);
julstep = C.sampling/1440; % minutes/(60*24)
julgrid = floor(julday(1)) : julstep : floor(julday(1))+1; julgrid = julgrid';
idxgrid = zeros(size(julgrid));
for ij = 1:length(julgrid);
    [minval,minindx] = min( abs(julday-julgrid(ij)) );
    if minval < julstep; idxgrid(ij) = minindx; end; 
end
% remove indices if empty
indx = find(idxgrid==0);
idxgrid(indx) = [];
% remove indices if double
indx = find(diff(idxgrid)==0);
idxgrid(indx) = [];
O.idxgrid = idxgrid;

% fill O
% TIME RESOLUTION IS REDUCED 1DVAR RET TIMES. MAYBE THIS IS NOT WANTED for ZH
O.angles_az = azi; % '0=North, 90=East, 180=South, 270=West'
O.angles_el = ele; % 'retrieval elevation angle'
O.channels = freq; % 'frequency of microwave channels [GHz]'
O.time = ( julday(idxgrid)-floor(julday(1)) ) * 24 * 3600; % seconds from midnight
O.y = tbs(:,idxgrid);  % 'brightness temperature [K]' - n_freq,time
O.ta = ta(idxgrid);  % 'air_temperature [K]'
O.pa = pa(idxgrid);  % 'air_pressure [Pa]'
O.hur = hur(idxgrid);  % 'relative_humidity [%/100]'
O.ybias = offset_tb(:,idxgrid); % n_freq,time
O.freq_shift = freq_shift;
if strcmp(C.instrument,'MP3000A')
   % This only applies to MP3000A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   O.wl_irp = wl_irp;          % 'infrared pyrometer central wavelength [micrometers]'
   O.tb_irp = tb_irp(idxgrid); % 'infrared pyrometer brightness temperature [K]'
   O.el_irp = el_irp(idxgrid); % 'infrared pyrometer elevation angle [degree]'
end

% building the R covariance matrix
% Dafault?
if C.Rdefault == 1; % 0/1 from_netcdf/default;
    switch C.instrument
        case 'HATPRO'
             % This corresponds to the one in Netcdf from Julich
             Rdiag = [1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    0.7000    0.5000    0.2000    0.2000    0.2000    0.2000]';
        case 'MP3000A'
             % This is adapted from the one above (Netcdf from Julich)
             Rdiag = [1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    0.7000    0.5000    0.2000    0.2000    0.2000    0.2000]';
    end
    O.Rcov = diag(Rdiag.^2); 
else
    O.Rcov = tb_cov + diag(tb_bias.^2); % it's missing the forward model uncertainty
end

% Diagonal?
if C.Rdiagonal == 1; % 0/1 full matrix/diagonal
   O.Rcov = diag(diag(O.Rcov));
end


return
