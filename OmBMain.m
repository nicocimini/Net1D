% O minus B monitoring: main program
%
% OmBMain reads the Config file, load observations and background, 
% line up Obs and Bkg, runs RTTOV-gb, store output.
% It loops on days from the selected period.

addpath('/Users/Nico/Matools/RTTOV/')

% Read configuration file (generates C structure)
Config_file = 'OmBConfig_000.m'; eval(Config_file(1:end-2)); Creset = C;

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

       % Load bkg data (includes definition of error covariance matrices)
       X = Net1DLoad_background(C);

       % Line up Obs and Bkg
       X = OmBLine_up_O_B(O,X);

       % Run RTTOV-gb
       X = OmBRun_RTTOV(X,C.instrument);
       
       % save output
       OmBSave_output(X,C);
       
    end
    
end

