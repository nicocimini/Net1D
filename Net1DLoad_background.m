% Network 1DVAR+RTTOV retrieval: Load background data
%
% Net1DLoad_background loads background arome data according to Config C 
% structure, i.e. for different station, day,...

function X = Net1DLoad_background(C);

% Initialization
X.time = []; % FixMe: null output
ij = C.bkgindx; % index of background profile (among the 9 extracted from Arome) - Olivier suggested #4, but since he uses python that means #5

% Searching for data files
bkgpath = [C.datapath 'arome/' C.station_id '/'];
bkghead = ['*ps_' C.station(1:3) '_arofr0*_l2_*'];
bkgdate = [num2str(C.day_one(1)) twodigstr(C.day_one(2)) twodigstr(C.day_one(3))];
listofiles = dir([bkgpath bkghead bkgdate '*']);

% If file do not exist, return
if isempty(listofiles); disp('No backgound files for this date:'); return; end;

% Just a check
if size(listofiles) ~= 3; disp('something is wrong with background files'); return; end;

% Load netcdf
% load ta and ancillary
bkgfile = [bkgpath listofiles(3).name];
lat = ncread(bkgfile,'lat',ij,1); % 'degree_north'                             - only reads n.ij (4 should be the center, according to Olivier, but Pauline thinks he was thinking in python, e.g. 4 is the fifth element)
lon = ncread(bkgfile,'lon',ij,1); % 'degree_east'                              
asl = ncread(bkgfile,'zsl',ij,1); % 'altitude above mean sea level [m]'        
ps = ncread(bkgfile,'ps',[ij 1],[1 Inf]);   % 'surface air_pressure [Pa]'      
ta = ncread(bkgfile,'ta',[ij 1 1],[1 Inf Inf]);        % 'air_temperature [K]' 
pa = ncread(bkgfile,'pressure',[ij 1 1],[1 Inf Inf]);  % 'air_pressure [Pa]'   
hg = ncread(bkgfile,'height',[ij 1 1],[1 Inf Inf]);    % height [m]            
% load hua
bkgfile = [bkgpath listofiles(1).name];
hua = ncread(bkgfile,'hua',[ij 1 1],[1 Inf Inf]); % 'absolute humidity [kg/m3]' 
% load hyd
bkgfile = [bkgpath listofiles(2).name];
clw = ncread(bkgfile,'clw',[ij 1 1],[1 Inf Inf]); % 'mass_fraction_of_clw_in_air [kg/kg??]'

% Select time
time = ncread(bkgfile,'time'); % forecast_time 'seconds since 1970-01-01 00:00:00 UTC'
time_ref = ncread(bkgfile,'forecast_reference_time'); % forecast_reference_time 'seconds since 1970-01-01 00:00:00 UTC'
hour = (time - time_ref)/3600;
indx3 = find(hour == 3); %% 3h-forecasts
%indx3 = [1; indx3]; % here I use the analysis to cover the first 3 hours (should be replaced by latest 3h forecast of previous day)
[timestr,julday,datetime] = computertime(time(indx3)); 

% Fill X
X.time = ( julday - floor(julday(1)) ) * 24 * 3600; % seconds from midnight;
X.lat = lat;
X.lon = lon;
X.asl = asl;
X.ps = ps(indx3);
% shall we interpolate on a fixed vertical grid before performing 1DVAR?
% Right now the interpolation is done before writing netcdf file
X.T = squeeze(ta(:,:,indx3));    
X.P = squeeze(pa(:,:,indx3));
X.Z = squeeze(hg(:,:,indx3));
X.Q = squeeze(hua(:,:,indx3));
X.LWC = squeeze(clw(:,:,indx3));
X.LWP=zeros(length(indx3),1);

% calculate liquid water content density
LWC_density=zeros(length(X.T(:,1)),1);
for i=1:length(indx3)
    LWC_density=mixingratio_to_density(X.LWC(:,i),X.T(:,i),X.P(:,i),X.Q(:,i));
    LWP_tempo=0;
    for lev=1:length(X.LWC(:,i))-1
        LWP_tempo=LWP_tempo+(LWC_density(lev)+LWC_density(lev+1))/2 * (X.Z(lev+1,i)-X.Z(lev,i));
    end
    X.LWP(i)=LWP_tempo;
end

% Flipping bottom-to-top to top-to-bottom to match 1DVAR input
X.T = flipud(X.T);
X.P = flipud(X.P);
X.Z = flipud(X.Z);
X.Q = flipud(X.Q);
X.LWC = flipud(X.LWC);
X.nlev = length(X.T(:,1));

% Add Bmatrix
switch C.Boptions
    
    case 0 % default
         C.Bfile = 'BTQ_OPER_FRANXL.txt';
    otherwise
         disp('Something is wrong with your B matrix'); return;        
end
X.B = load([C.Bpath C.Bfile]);

end
                      