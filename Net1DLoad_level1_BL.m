% Network 1DVAR+RTTOV retrieval: Load level1 BL data
%
% Net1DLoad_level1_BL loads BL level 1 files according to Config C structure, i.e. for
% different station, instrument, day,...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HATPRO and MP3000A BL files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [O,C] = Net1DLoad_level1_BL(C,lv1file);

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
try; O.GA.Processing_date = ncreadatt(lv1file,'/','Processing_date'); catch; O.GA.Processing_date = ' '; end;
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
        % Nico 2017/05/25
        % moved this within the "if strcmp(C.instrument,'MP3000A')" statement
        % i.e. after the elevation angle resampling
        %badindx = find(flag);
    case 'sir'
        badindx = find(not(isnan(flag)));
end


% load data
time = ncread(lv1file,'time'); % 'seconds since 1970-01-01 00:00:00 UTC'
azi = ncread(lv1file,'azi'); % '0=North, 90=East, 180=South, 270=West'
ele = ncread(lv1file,'ele'); % 'retrieval elevation angle'
freq = ncread(lv1file,'freq_sb'); % 'frequency of microwave channels [GHz]'
tbs = ncread(lv1file,'tb');  % 'brightness temperature [K]' - n_freq,n_angle,time
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


% This only applies to MP3000A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(C.instrument,'MP3000A')
    elangs = unique(ele);
    indx15 = find(ele==elangs(1));
    indx90 = find(ele==elangs(2));
    if size(elangs) ~= 2; disp('something is wrong with MP3000A scan'); return; end;
    if size(indx15) ~= size(indx90); disp('something is wrong with MP3000A scan'); return; end;
    time = time(indx15);
    flag = flag(indx15) + flag(indx90); % accounting for both angles (0: none is flagged, 1: one is flagged, 2: both are flagged)
    badindx = find(flag);
    azi = unique(azi);
    if length(azi)==1 & length(elangs)>1; azi = azi * ones(size(elangs)); end;
    ele = flipud(elangs); % put 90? degrees first
    tbs2(:,1,:) = tbs(:,indx90);
    tbs2(:,2,:) = tbs(:,indx15);
    tbs = tbs2; clear tbs2;
    ta = ta(indx15);
    pa = pa(indx15);
    hur = hur(indx15);
    offset_tb2(:,1,:) = offset_tb(:,indx90);
    offset_tb2(:,2,:) = offset_tb(:,indx15);
    offset_tb = offset_tb2; clear offset_tb2;
    % load extra fields
    wl_irp = ncread(lv1file,'wl_irp'); % 'infrared pyrometer central wavelength [micrometers]'
    tb_irp = ncread(lv1file,'tb_irp'); % 'infrared pyrometer brightness temperature [K]'
    el_irp = ncread(lv1file,'ele_irp'); % 'infrared pyrometer elevation angle [degree]'
    tb_irp = tb_irp(indx15)';
    el_irp = el_irp(indx15);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Other instrument specific fix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To be removed one day we make the provider fix these %%%%%%%%%%%%%%%%%%%%
switch C.station_id
    case 'ces'
          tbs = permute(tbs,[2 1 3]);             % original files switch the nchn and nang dimensions
          offset_tb = permute(offset_tb,[2 1 3]); % original files switch the nchn and nang dimensions
    case 'joy'

    case 'lac'

    case 'pay'

    case 'rao'

    case 'sir'
          tbs = tbs(:,1:6,:); % discarding the observations at 0 elevation 
          ele = ele(1:6);     % discarding the observations at 0 elevation
          ele = ele - 9*1e4;  % don't know the reason of this constant;
          %ele(1:6) = ele(1:6) - 9*1e4; % don't know the reason...
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discard flagged data
time(badindx) = []; 
tbs(:,:,badindx) = []; 
ta(badindx) = [];
pa(badindx) = [];
hur(badindx) = [];
offset_tb(:,:,badindx) = []; 
if strcmp(C.instrument,'MP3000A')
    tb_irp(badindx) = [];
    el_irp(badindx) = [];
end

% Select time
%time_bnds = ncread(lv1file,'time_bnds');
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

% Here it computes the cloud flag based on std(Tb31) at the reduced resolution
julstd = C.std31wndw/1440; % minutes/(60*24)
switch C.instrument
    case 'HATPRO'; 
        chn31 = 7;
        chn23 = 2;
    case 'MP3000A'; 
        chn31 = 5; % ChECK IT! 
        chn23= 2;
end
ngrid = length(idxgrid);
std31 = ones(ngrid,1);
for in = 1:ngrid    
    idx = idxgrid(in);
    idxwndw = find( abs (julday - julday(idx)) < julstd );    
    std31(in) = std( tbs(chn31,1,idxwndw) ); % std(Tb31(90))   
end
cld31 = not ( not ( sign( std31-C.std31value ) + 1 ) ); % cld flag based on std

% fill O
O.angles_az = azi; % '0=North, 90=East, 180=South, 270=West'
O.angles_el = ele; % 'retrieval elevation angle'
O.channels = freq; % 'frequency of microwave channels [GHz]'
O.time = ( julday(idxgrid)-floor(julday(1)) ) * 24 * 3600; % seconds from midnight
O.y = tbs(:,:,idxgrid);  % 'brightness temperature [K]' - n_freq,n_angle,time
O.ta = ta(idxgrid);  % 'air_temperature [K]'
O.pa = pa(idxgrid);  % 'air_pressure [Pa]'
O.hur = hur(idxgrid);  % 'relative_humidity [%/100]'
%O.ybias = offset_tb(:,:,idxgrid); 
%O.yoffs = offset_tb(:,:,idxgrid); % this results from offset_tb (see description above) - I don't think we need this, as it comes from site dependent O-B
O.ybias = tb_bias;                % this results from tb_bias (see description above) - this shall be used for the systematic uncertainty of Tb
O.freq_shift = freq_shift;
if strcmp(C.instrument,'MP3000A')
   % This only applies to MP3000A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   O.wl_irp = wl_irp;          % 'infrared pyrometer central wavelength [micrometers]'
   O.tb_irp = tb_irp(idxgrid); % 'infrared pyrometer brightness temperature [K]'
   O.el_irp = el_irp(idxgrid); % 'infrared pyrometer elevation angle [degree]'
end
O.std31 = std31; % std of Tb31 (over C.std31wndw minutes)
O.cld31 = cld31; % cld flag based on std31

% Pauline compute LWP value from dual channel algorithm from Nico's
% function
LWP_values=zeros(ngrid,1);
for i=1:ngrid
    [V,L]=Net1DComp_LWP(tbs(chn23,1,idxgrid(i)),tbs(chn31,1,idxgrid(i)));
    LWP_values(i)=L;
end
O.LWP=LWP_values;

% Overwriting O.ybias if set to default.
% From HDCP2 netcedf file comment: 
% This variable is an estimate of the one-standard-deviation calibration error 
% to be expected from an absolute system calibration, i.e. the likely systematic 
% error of brightness temperature. As a reference see Maschwitz et al. 2013, AMT (Tab. 5). 
% However, these numbers differ from instrument to instrument and should be adapted accordingly. 
% Values only valid for elevation angles larger than 20deg.
if C.Usys_default
    switch C.instrument
        case 'HATPRO'
             % This corresponds to the one in Netcdf from Julich
             % O.ybias = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.7 0.5 0.2 0.2 0.2 0.2]';
             % This corresponds to Maschwitz, 2013 (Table 5). I take LN2
             % values because all considered HATPROs do not use tip cal.
             O.ybias = [1.2 1.6 1.0 0.9 1.1 1.0 1.0 1.0 0.7 0.5 0.2 0.3 0.3 0.2]';
        case 'MP3000A'
             % tb_bias is Nan in rao lv1 files
             % This is adapted from the Netcdf from Julich (see above)
             % O.ybias = [1.0 1.0 1.0 1.0 1.0 1.0 0.7 0.5 0.2 0.2 0.2 0.2]';
             % This is from Hewison 2006 (Table 3-2). I used the uncertainty 
             % of known biases (assuming they are removed by the instrument software): 
             % channels 1-5: TIP (Table 3-2, 6th column); channels 6-12: LN2 (Table 3-2, 4th column) 
             O.ybias = [0.52 0.40 0.40 0.29 0.21 1.06 0.89 0.38 0.24 0.22 0.22 0.22]';
    end    
end

% Building the R covariance matrix
% Following Hewison 2006, R = E + F + M, where: 
% E: measurement noise, F: modelling error, M: representativeness error
% Default E, F, and M are taken from Hewison 2006. 
% For site-specific values: E = tb_cov; while F and M should be estimated dinamically.
% These correspond to Hewison 2006 (Table 4-1), evaluated on a MP3000
Esigma = [0.17 0.12 0.11 0.13 0.21 0.18 0.15 0.17 0.18 0.19 0.54 0.18]';
Fsigma = [0.83 0.84 0.82 0.67 0.61 1.10 0.88 0.35 0.06 0.05 0.05 0.06]';
Msigma = [0.65 0.67 0.69 0.78 1.00 1.70 1.35 0.32 0.10 0.10 0.40 0.11]';
if C.instrument == 'HATPRO'
    % This corresponds to the one in Netcdf from Julich
    % Rsigma = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.7 0.5 0.2 0.2 0.2 0.2]';
    % This is derived from sqrt(diag(tb_cov)) from joy
    Esigma = [0.12 0.11 0.10 0.09 0.10 0.10 0.10 0.42 0.46 0.26 0.26 0.18 0.18 0.59]';
    % These are interpolated on HATPRO chan from MP3000 chan taken from Hewison 2006 (Table 4-1)
    frq_HATPRO = [22.24 23.04 23.84 25.44 26.24 27.84 31.40 51.26 52.28 53.86 54.94 56.66 57.30 58.00]';
    frq_MP3000 = [22.235 23.035 23.835 26.235 30.00 51.25 52.28 53.85 54.94 56.66 57.29 58.80]';
    Fsigma = interp1(frq_MP3000,Fsigma,frq_HATPRO);
    Msigma = interp1(frq_MP3000,Msigma,frq_HATPRO);
end
E = diag(Esigma.^2);
F = diag(Fsigma.^2);
M = diag(Msigma.^2);
if ~C.Rdefault % 0/1 for netcdf/default;
    E = tb_cov; % [K^2]
end
O.Rcov = E + F + M; 

% Forcing R to be diagonal, if required
if C.Rdiagonal % 0/1 for full-matrix/diagonal
   O.Rcov = diag(diag(O.Rcov));
end


end
