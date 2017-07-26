% Network 1DVAR+RTTOV retrieval (Net1D): configuration file
%
% Net1DConfig fills the C structure, which is used for all the computations.
% The C structure shall be kept with the 1DVAR+RTTOV output.
%
% Use this file to change all the settings (paths, station, period, bias
% correction, ...). 
% Change the progressive number (Net1DConfig_NNN.m) to keep the used Config file

% Make sure C is empty
clear C

% Paths
%C.datapath = '/Users/Nico/PROGETTI/IMAA/GAIACLIM/DATA/VO_MWR/TOPROF/'; % 2015 dataset
%C.datapath = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/'; % 2014 dataset

%directory with the input data
C.datapath='/home/martinet/1DVAR/NET1D/DATA/TOPROF_DATA/TOPROF/';
%C.ODVARpath = '/home/martinet/1DVAR/NWP_SAF/1DVAR_MWR/1DVar_K_complete_clean_v1.2-debugQ/';
% directory with the 1DVAR executable
C.ODVARpath= '/home/martinet/1DVAR/NWP_SAF/1DVAR_MWR/SAF1DVARgb/SAF1DVARgb_GIT/NWPSAF_1DVar_1.1/1DVar/';
C.ODVARpath_retrieval='/home/martinet/1DVAR/NWP_SAF/1DVAR_MWR/SAF1DVARgb/SAF1DVARgb_GIT/NWPSAF_1DVar_1.1/1DVar/Output_MWR/';
%directory where to save the output .mat and netcdf 1DVAR files
C.ODVARpath_output = [C.datapath '/1DVAR/Config_001/']; % main 1DVAR folder

%path to RTTOV coefficient files
C.rttocoefpath='/home/martinet/rttov/rttov_visee_sol/coeffiles/';

% Station, instrument, time, 
C.station = 'payerne'; C.instrument = 'HATPRO'; C.channum = 14;
%C.station = 'rao'; C.instrument = 'MP3000A'; C.channum = 12;
C.day_start = [2014 1 1]; % YYYY MM dd
C.day_end =   [2014 1 1]; % YYYY MM dd
C.sampling = 60; % minutes


% Instrument config
C.elev_angles_how_many = 1;   % >=1
C.elev_angles_degrees = [90 42]; % size of C.elev_angles_how_many
C.elev_angles_channum{1} = [1:C.channum]; % size of C.elev_angles_how_many
C.elev_angles_channum{2} = [11:14];

% Bias correction options
C.biascorrection = 1; % 0/1
if C.biascorrection
   C.biascorr_path = '/home/martinet/1DVAR/NET1D/DATA/TOPROF_DATA/OB_bias/'; % where to get bias to be applied
   C.biascorr_file = ['bias_' C.station(1:3) '_000.mat']; % add another file for updating bias correction
end

% Control variables
C.retrieve_T = [1 1 60]; % [T/F "1st-lev within bkg profile" "how many lev from 1st-lev"]
C.retrieve_Q = [1 1 60]; % 
C.retrieve_LWP = [1]; % 
%C.retrieve_LWC = []; % 


% Covariance error matrix
C.Rdefault = 1; % 0/1 from_netcdf/default;
C.Rdiagonal = 1; % 0/1 full matrix/diagonal
C.Roptions = [0/1 C.Rdiagonal]; % options for R:
C.R = []; % Default R; number of channels
C.Rcov = zeros(); % Default Rcov;
C.Boptions = 0; % options for B: 0/1/2.. default/something else;
C.Bpath = '/home/martinet/1DVAR/NET1D/DATA/TOPROF_DATA/Bmatrix/'; % path B (default or others);


% Sky conditions
C.Skyconditions = 'all'; % process clear sky, cloudy sky, or other options
C.force_bgk_clear = 0;   % force the bacground LWC to zero

% Output format and other characteristics
C.Output_format = 'matlab'; % either 'matlab' or 'netcdf' (for future development)


% Minimization options
C.Minimisation_Method = 'ML'; % either 'ML' for Marquard-Levenberg or 'N' for Newtonian
C.MaxIterations = 10; % Maximum number of iterations in the external loop:
C.RTTOV_version = 'RTTOV11'; % version of RTTOV used for the retrieval: for now only 'RTTOV11' is available
C.GeneralMode_1DVAR = 40; % 0='Operational', 10='Production', 20='Diagnostic', 30='Debug', 40='Verbose'.
C.MissingValue = -999;

% Station info
C.station_id = C.station(1:3);
switch C.station_id
    case 'ces'
          C.lat = 51.97; C.lon = 4.93; C.asl = -0.7;
    case 'joy'
          C.lat = 50.91; C.lon = 6.41; C.asl = 111;
    case 'lac'
          C.lat = 51.35; C.lon = 12.43; C.asl = 125;
    case 'pay'
          C.lat = 46.82; C.lon = 6.95; C.asl = 491;
    case 'rao'
          C.lat = 52.21; C.lon = 14.12; C.asl = 125;
    case 'sir'
          C.lat = 48.80; C.lon = 2.36; C.asl = 156;
end

% Section for Netcdf output file attributes
% % Units
% C.nc.units.time = 'seconds since 1970-01-01 00:00:00 UTC';
% C.nc.units.freq_sb = 'GHz';
% C.nc.units.lat = 'degree_north';
% C.nc.units.lon = 'degree_east';
% C.nc.units.zsl = 'm';
% % Standard name
% C.nc.standard_name.time = 'time';
% C.nc.standard_name.freq_sb = 'sensor_band_central_radiation_frequency';
% C.nc.standard_name.lat = 'latitude';
% C.nc.standard_name.lon = 'longitude';
% C.nc.standard_name.zsl = 'altitude';
% % Long name
% C.nc.long_name.freq_sb = 'frequency of microwave channels';
% C.nc.long_name.zsl = 'altitude above mean sea level';
% % Bunds
% C.nc.bounds.time = 'time_bnds';
% Global attributes
C.nc.glatt.Institution      = 'CNR-IMAA';
C.nc.glatt.Contact_person   = 'Reprocessed data: Nico Cimini, CNR-IMAA (domenico.cimini@imaa.cnr.it). Original data: ';
C.nc.glatt.History          = 'Data reprocessed with Net1D v1.0. Original ';
C.nc.glatt.Processing_date  = datestr(now);
C.nc.glatt.Author           = 'Reprocessed data: Nico Cimini, CNR-IMAA (domenico.cimini@imaa.cnr.it). Original data: ';
C.nc.glatt.Comments         = 'Reprocessed data for the GAIA-CLIM Virtual Observatory from original data. ';


% Fixed Z grid (top to bottom)
%C.FixZgrid = [46000 33000 27000 23000:-2000:17000 16000:-1000:10000 9700:-600:7900 7400:-500:6400 6100:-300:3900 3700:-200:1900 1750:-150:1000 900:-100:400 350:-60:170 130:-40:50 30:-20:10]';
C.FixZgrid = [10000:-500:5000 4750:-250:2000 1900:-100:1000 950:-50:200 175:-25:50 40:-10:10]';
