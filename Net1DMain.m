% Network 1DVAR+RTTOV retrieval: main program
%
% Net1DMain reads the Config file, load observations and background, apply
% bias correction (if required), runs 1DVAR, store output.
% It loops on days from the selected period.
%
% Net1D development notes:
% https://pad.okfn.org/p/Net1D_development

% Read configuration file (generates C structure)
Config_file = 'Net1DConfig_000.m'; eval(Config_file(1:end-2)); Creset = C;

% Save the configuration file in the output directory
if ~exist(C.Configoutpath,'dir'); mkdir(C.Configoutpath); end
command = ['cp ' Config_file ' ' C.Configoutpath]; system(command);


% Loop on days
% C.day_one = C.day_start; % this has to go from start to end
%dayindx = [datenum(C.day_start) : datenum(C.day_end)];
for id = C.dayindx' % must be row vector

    % Set date
    C = Creset; % reset C
    day_one = datevec(id);
    C.day_one = day_one(1:3);
    disp(['Processing ' C.station ' ' C.instrument ' for ' datestr(id)]);

    % Load lv1 data (includes definition of error covariance matrices)
    [O,C] = Net1DLoad_level1(C); 

    % Throw message
    disp([O.content ' file found']);

    % Continue only if lv1 file was found
    if ~strcmp(O.content,'No lv1'); 

       % Apply bias correction
       O = Net1DBias_correction(C,O);

       % Load bkg data (includes definition of error covariance matrices)
       X = Net1DLoad_background(C);

       % Run 1DVAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Clean 1DVAR directory and link to RTTOV coeff file
       Net1DClean_1DVAR(C)
       % Write input 1DVAR
       X = Net1DWrite_1DVAR_MAIN(O,X,C);
       % Run 1DVAR
       thisdir = pwd;
       cd([C.ODVARpath '/WorkDir/']);
       command = ['./Run_1DVar_test_MWR.ksh ' C.instrument];
       [status,result] = system(command);
       cd(thisdir);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       % Retrieve output and compute error budget
       R = Net1DLoad_1DVARout(C,X);
       AK = Net1DLoad_1DVARout_AK(C,X);
       J = Net1DLoad_1DVARout_Jacobians(C,X);
       A = Net1DLoad_1DVARout_A(C,X);
       E = Net1DComp_errorbudget(C,O.Rcov,X.B,AK,J);
    
       % Save output (write Matlab, ASCII, and/or Netcdf files)
       Net1DSave_1DVARout(C,O,X,R,A,E);
    
       % Plotting (just for checking)
       Net1DPlot_1DVARout(C,O,X,R,E,A,AK,J);
    
    end
    
end
