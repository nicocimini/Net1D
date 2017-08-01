% Network 1DVAR+RTTOV retrieval: main program
%
% Net1DMain reads the Config file, load observations and background, apply
% bias correction (if required), runs 1DVAR, store output.
% It loops on days from the selected period.
%
% Net1D development notes:
% https://pad.okfn.org/p/Net1D_development

% Read configuration file (generates C structure)
Config_file = 'Net1DConfig_000';
eval(Config_file); C.Config_file = [Config_file '.m'];

%%%%save the configuration file in the output directory
out_path = [C.ODVARpath_output 'Config/' C.station_id '/'];
if ~exist(out_path,'dir'); mkdir(out_path); end
command = ['cp ' C.Config_file ' ' out_path]; system(command);


% FixMe: here should go the loop on days (take first day for now)
C.day_one = C.day_start; % this has to go from start to end


% Load data (includes definition of error covariance matrices)
[O,C] = Net1DLoad_level1(C);
X = Net1DLoad_background(C);

% Apply bias correction
O = Net1DBias_correction(C,O);


% Run 1DVAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clean 1DVAR directory and link to RTTOV coeff file
Net1DClean_1DVAR(C)

% Write input 1DVAR
X = Net1DWrite_1DVAR_MAIN(O,X,C);

% run 1DVAR
% either a matlab or a shell script (simbolic link or copy files) -> Pauline
thisdir = pwd;
cd([C.ODVARpath '/WorkDir/']);
command = ['./Run_1DVar_test_MWR.ksh ' C.instrument];
[status,result] = system(command);
cd(thisdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% retrieve output and compute error budget
R = Net1DLoad_1DVARout(C,X);
AK = Net1DLoad_1DVARout_AK(C,X);
J = Net1DLoad_1DVARout_Jacobians(C,X);
A = Net1DLoad_1DVARout_A(C,X);
E = Net1DComp_errorbudget(C,O.Rcov,X.B,AK,J);

% save output
% write output in either Matlab or Netcdf format
Net1DSave_1DVARout(C,O,X,R,A,E); 

% Plotting (just for checking)
Net1DPlot_1DVARout(C,O,X,R,E,A,AK,J);
