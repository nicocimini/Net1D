% O minus B monitoring: Run RTTOV-gb
%
% It writes the input, runs RTTOV-gb, and retrieves the output
%
% This is highly inefficient, as it writes input, call RTTOV, and reads
% output for each profile and each elev angle. To be improved in the future, i.e.
% to be replaced with one call for all profiles; need to modify EXAMPLE_FWD.F90

function X = OmBRun_RTTOV(X,instrument)

% Set paths
path_RTTOV = '/Users/Nico/SFTWR/RTTOV/RTTOVgb_v1.3/rttov112_GB_K_v1.3/rttov_test/nico_matlab_call/';

% Set dimensions
ntime = length(X.time);
neang = length(X.angles_el);
nchan = length(X.channels);
TbBkg = zeros(size(X.TbObs));

% Start loop on profiles 
% This is highly inefficient, to be replaced with one call for all profiles (need modify EXAMPLE_FWD.F90)
for it = 1:ntime
    for ia = 1:neang
        
        % Write RTTOV-gb input
        Write_RTTOVgb_input(path_RTTOV,X.P(:,it),X.T(:,it),X.Q(:,it),90-X.angles_el(ia));
        
        % Run RTTOV-gb
%        system([path_RTTOV 'run_matlab_call_' instrument '.sh']);
        thisdir = pwd;
        cd(path_RTTOV)
        system(['./run_matlab_call_' instrument '.sh']);
        cd(thisdir)
        
        % Retrieve RTTOV-gb output
        TbBkg(:,ia,it) = Load_RTTOVgb_output([path_RTTOV 'output_example_fwd.dat'],nchan);

    end
end

% Add info to X
X.instrument = instrument;
X.TbBkg = TbBkg;


return