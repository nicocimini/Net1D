% Network 1DVAR+RTTOV retrieval: Save 1DVAR output data in netcdf
%
% Net1DSave_netcdflv2 saves output 1DVAR data according to Config C structure
% The output file is a CF-compliant netcdf file, as checked through: 
% http://puma.nerc.ac.uk/cgi-bin/cf-checker.pl

% lv2file = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/level2/joy/sups_joy_mwrBL00_l2_ta_v01_20140101005346.nc'
% sups_joy_mwr00_l2_clwvi_i00_20140101000004.nc	
% sups_joy_mwr00_l2_hua_i00_20140101000004.nc	
% sups_joy_mwrBL00_l2_ta_v01_20140101005346.nc
% sups_joy_mwr00_l2_prw_i00_20140101000004.nc	
% ncdisp(lv2file);

function Net1DSave_netcdflv2(l1outputfile,C,O,X,R,A,E);


% Write ta file
% lv2file = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/level2/joy/sups_joy_mwrBL00_l2_ta_v01_20140101005346.nc'
% ncdisp(lv2file);
if C.retrieve_T(1)
   Net1DSave_netcdflv2_ta(l1outputfile,C,O,X,R,A,E);
end

% Write hua file
% lv2file = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/level2/joy/sups_joy_mwr00_l2_hua_i00_20140101000004.nc'
% ncdisp(lv2file);
if C.retrieve_Q(1)
   Net1DSave_netcdflv2_hua(l1outputfile,C,O,X,R,A,E);
   Net1DSave_netcdflv2_prw(l1outputfile,C,O,R,E);
end

% Write clwvi file
% lv2file = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/level2/joy/sups_joy_mwr00_l2_clwvi_i00_20140101000004.nc'
% ncdisp(lv2file);


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub-functions to write different NetCDF file types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Net1DSave_netcdflv2_prw %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Net1DSave_netcdflv2_prw(l1outputfile,C,O,R,E);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change output name (to do for each variable/file, e.g. hua, ipwv,...)
ncoutputfile = strrep(l1outputfile,'tb','prw');
if exist(ncoutputfile); delete(ncoutputfile); end; % delete lv2 netcdf file if already exists

% Select variables to write
% Variables to be written in the lev2 netcdf file
vars = strvcat('time','time_bnds','lat','lon','zsl','azi','ele','ele_ret','prw','prw_err','prw_sys','flag','prw_offset','prw_off_zenith','prw_off_zenith_offset');

% Converting seconds from midnight to seconds from 1/1/1970
time_bnds = C.sampling;
time_offset = ( datenum(C.day_one(1),C.day_one(2),C.day_one(3)) - datenum(1970,1,1) ) * 24 * 3600;

% Interpolate on a fixed height grid
nret = length(R);
time = zeros(nret,1); 
prw     = C.MissingValue * ones(nret,1); % set to missing value
prw_err = C.MissingValue * ones(nret,1); % set to missing value
prw_sys = C.MissingValue * ones(nret,1); % set to missing value
for ir = 1:nret
    time(ir) = R(ir).time + time_offset;
    if R(ir).nite <= C.MaxIterations % convergence was reached
       prw(ir)     = R(ir).Ret_IWV_kgm2;  % [kg/m2]==[mm] - TWVC value
       prw_err(ir) = E(ir).IWV_rnd;       % [kg/m2]==[mm] - TWVC rand uncertainty (based on ObsErr only)
       prw_sys(ir) = E(ir).IWV_sys;       % [kg/m2]==[mm] - TWVC syst uncertainty (summed linearly)
    end
end
flag                  = zeros(nret,1);                 % I leave these fields unused unless they are needed in the future
prw_offset            = C.MissingValue * ones(nret,1); % "
prw_off_zenith        = C.MissingValue * ones(nret,1); % "
prw_off_zenith_offset = C.MissingValue * ones(nret,1); % "

% Create variables  
nccreate(ncoutputfile,'time','Dimensions',{'time' length(time)},'Datatype','double','Format','netcdf4');
nccreate(ncoutputfile,'time_bnds','Dimensions',{'nv' 2 'time'},'Datatype','double');
nccreate(ncoutputfile,'lat','Datatype','single');
nccreate(ncoutputfile,'lon','Datatype','single');
nccreate(ncoutputfile,'zsl','Datatype','single');
nccreate(ncoutputfile,'azi','Dimensions',{'n_angle' length(O.angles_az)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'ele','Dimensions',{'n_angle' length(O.angles_el)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'ele_ret','Dimensions',{'time'},'Datatype','single');                       % this may have to change for RAO!
nccreate(ncoutputfile,'prw','Dimensions',{'time'},'Datatype','single');  
nccreate(ncoutputfile,'prw_err','Dimensions',{'time'},'Datatype','single'); % different than original lv2 (pwr_err has n_ret dimention)
nccreate(ncoutputfile,'prw_sys','Dimensions',{'time'},'Datatype','single'); % new wrt to original lv2. pwr_sys has time dimension
nccreate(ncoutputfile,'flag','Dimensions',{'time'},'Datatype','single');                  % unused
nccreate(ncoutputfile,'prw_offset','Dimensions',{'time'},'Datatype','single');            % unused
nccreate(ncoutputfile,'prw_off_zenith','Dimensions',{'time'},'Datatype','single');        % unused
nccreate(ncoutputfile,'prw_off_zenith_offset','Dimensions',{'time'},'Datatype','single'); % unused

% Write data associated to variables
ncwrite(ncoutputfile,'time',time);
ncwrite(ncoutputfile,'time_bnds',[time-time_bnds time+time_bnds]');
ncwrite(ncoutputfile,'lat',C.lat);
ncwrite(ncoutputfile,'lon',C.lon);
ncwrite(ncoutputfile,'zsl',C.asl);
ncwrite(ncoutputfile,'azi',O.angles_az);
ncwrite(ncoutputfile,'ele',O.angles_el);
ncwrite(ncoutputfile,'ele_ret',90);             % We always consider retrievals at zenith, right?
ncwrite(ncoutputfile,'prw',prw); 
ncwrite(ncoutputfile,'prw_err',prw_err);
ncwrite(ncoutputfile,'prw_sys',prw_sys);
ncwrite(ncoutputfile,'flag',flag);                                  % unused
ncwrite(ncoutputfile,'prw_offset',prw_offset);                      % unused
ncwrite(ncoutputfile,'prw_off_zenith',prw_off_zenith);              % unused
ncwrite(ncoutputfile,'prw_off_zenith_offset',prw_off_zenith_offset);% unused

% Set variable attributes
Variables = Net1DVariables_lv2;
nvar = length(Variables);
for iv = 1:nvar
    varName = Variables(iv).Name;  
    if ismember(varName,vars,'rows')        
              ncwriteallatt(ncoutputfile,Variables(iv));
    end    
end

% Set global attributes
ncwriteatt(ncoutputfile,'/','Title','Microwave radiometer retrieved precipitable water (prw), also called total water vapor content (TWVC)');
ncwriteatt(ncoutputfile,'/','Institution',[O.GA.Institution]);
ncwriteatt(ncoutputfile,'/','Contact_person',[C.nc.glatt.Contact_person O.GA.Contact_person]);
ncwriteatt(ncoutputfile,'/','Source',[O.GA.Source]);
ncwriteatt(ncoutputfile,'/','History',[C.nc.glatt.History]);
ncwriteatt(ncoutputfile,'/','Dependencies',l1outputfile);
ncwriteatt(ncoutputfile,'/','Conventions',[O.GA.Conventions]);
ncwriteatt(ncoutputfile,'/','Processing_date',[C.nc.glatt.Processing_date]);
ncwriteatt(ncoutputfile,'/','Author',[C.nc.glatt.Author O.GA.Author]);
ncwriteatt(ncoutputfile,'/','Comments',[C.nc.glatt.Comments O.GA.Comments]);
ncwriteatt(ncoutputfile,'/','License',[O.GA.License]);
ncwriteatt(ncoutputfile,'/','Measurement_site',[O.GA.Measurement_site]);

return

% Net1DSave_netcdflv2_hua %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Net1DSave_netcdflv2_hua(l1outputfile,C,O,X,R,A,E);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change output name (to do for each variable/file, e.g. hua, ipwv,...)
ncoutputfile = strrep(l1outputfile,'tb','hua');
if exist(ncoutputfile); delete(ncoutputfile); end; % delete lv2 netcdf file if already exists

% Select variables to write
% Variables to be written in the lev2 netcdf file
vars = strvcat('time','time_bnds','lat','lon','zsl','azi','ele','ele_ret','height','hua','hua_offset','hua_err','hua_sys');

% Converting seconds from midnight to seconds from 1/1/1970
% NB: Be careful here: does O.time and R matches? How about non convergence episodes? 
% If R(i).nite >= C.MaxIterations means no convergence. Should treat this!
%It should be OK now, as time is filled within the loop for ir = 1:nret below
time_bnds = C.sampling;
time_offset = ( datenum(C.day_one(1),C.day_one(2),C.day_one(3)) - datenum(1970,1,1) ) * 24 * 3600;

% Interpolate on a fixed height grid
nret = length(R);
nlev = length(C.FixZgrid);
time = zeros(nret,1); 
hua     = C.MissingValue * ones(nlev,nret); % set to missing value
hua_err = C.MissingValue * ones(nlev,nret); % set to missing value
hua_sys = C.MissingValue * ones(nlev,nret); % set to missing value
for ir = 1:nret
    time(ir) = R(ir).time + time_offset;
    if R(ir).nite <= C.MaxIterations % convergence was reached
       ibkg = X.ObsMatch(R(ir).nobs_prf);
       if ( X.Z(1,ibkg) < C.FixZgrid(1) ) | ( X.Z(end,ibkg) > C.FixZgrid(end) ); 
          disp('Interpolation may produce NaN at the vertical boudaries'); return; 
       end
       hua(:,ir)     = interp1( X.Z(:,ibkg) , R(ir).Ret_H_kgm3 , C.FixZgrid);
       hua_err(:,ir) = interp1( X.Z(:,ibkg) , A(ir).HTotErr    , C.FixZgrid); 
       hua_sys(:,ir) = interp1( X.Z(:,ibkg) , E(ir).HSysUnc    , C.FixZgrid); 
    end
end
hua_offset = zeros(size(hua_err(:,1))); % I leave this field zero unless it's needed in the future


% Create variables  
nccreate(ncoutputfile,'time','Dimensions',{'time' length(time)},'Datatype','double','Format','netcdf4');
nccreate(ncoutputfile,'time_bnds','Dimensions',{'nv' 2 'time'},'Datatype','double');
nccreate(ncoutputfile,'lat','Datatype','single');
nccreate(ncoutputfile,'lon','Datatype','single');
nccreate(ncoutputfile,'zsl','Datatype','single');
nccreate(ncoutputfile,'azi','Dimensions',{'n_angle' length(O.angles_az)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'ele','Dimensions',{'n_angle' length(O.angles_el)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'ele_ret','Dimensions',{'time'},'Datatype','single');                       % this may have to change for RAO!
nccreate(ncoutputfile,'height','Dimensions',{'height' nlev},'Datatype','single');
nccreate(ncoutputfile,'hua','Dimensions',{'height' 'time'},'Datatype','single');  
nccreate(ncoutputfile,'hua_offset','Dimensions',{'height'},'Datatype','single');     % different than original lv2 (hua_offset has height and time dimension)
nccreate(ncoutputfile,'hua_sys','Dimensions',{'height' 'time'},'Datatype','single'); % new wrt to original lv2. hua_sys has height and time dimension
nccreate(ncoutputfile,'hua_err','Dimensions',{'height' 'time'},'Datatype','single'); % different than original lv2 (hua_err has only height dimension)

% Write data associated to variables
ncwrite(ncoutputfile,'time',time);
ncwrite(ncoutputfile,'time_bnds',[time-time_bnds time+time_bnds]');
ncwrite(ncoutputfile,'lat',C.lat);
ncwrite(ncoutputfile,'lon',C.lon);
ncwrite(ncoutputfile,'zsl',C.asl);
ncwrite(ncoutputfile,'azi',O.angles_az);
ncwrite(ncoutputfile,'ele',O.angles_el);
ncwrite(ncoutputfile,'ele_ret',90);             % We always consider retrievals at zenith, right?
ncwrite(ncoutputfile,'height',C.FixZgrid);
ncwrite(ncoutputfile,'hua',hua); 
ncwrite(ncoutputfile,'hua_offset',hua_offset); 
ncwrite(ncoutputfile,'hua_err',hua_err);
ncwrite(ncoutputfile,'hua_sys',hua_sys);

% Set variable attributes
Variables = Net1DVariables_lv2;
nvar = length(Variables);
for iv = 1:nvar
    varName = Variables(iv).Name;  
    if ismember(varName,vars,'rows')        
              ncwriteallatt(ncoutputfile,Variables(iv));
    end    
end

% Set global attributes
ncwriteatt(ncoutputfile,'/','Title','Microwave radiometer retrieved humidity profile');
ncwriteatt(ncoutputfile,'/','Institution',[O.GA.Institution]);
ncwriteatt(ncoutputfile,'/','Contact_person',[C.nc.glatt.Contact_person O.GA.Contact_person]);
ncwriteatt(ncoutputfile,'/','Source',[O.GA.Source]);
ncwriteatt(ncoutputfile,'/','History',[C.nc.glatt.History]);
ncwriteatt(ncoutputfile,'/','Dependencies',l1outputfile);
ncwriteatt(ncoutputfile,'/','Conventions',[O.GA.Conventions]);
ncwriteatt(ncoutputfile,'/','Processing_date',[C.nc.glatt.Processing_date]);
ncwriteatt(ncoutputfile,'/','Author',[C.nc.glatt.Author O.GA.Author]);
ncwriteatt(ncoutputfile,'/','Comments',[C.nc.glatt.Comments O.GA.Comments]);
ncwriteatt(ncoutputfile,'/','License',[O.GA.License]);
ncwriteatt(ncoutputfile,'/','Measurement_site',[O.GA.Measurement_site]);

return

% Net1DSave_netcdflv2_ta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Net1DSave_netcdflv2_ta(l1outputfile,C,O,X,R,A,E);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change output name (to do for each variable/file, e.g. hua, ipwv,...)
ncoutputfile = strrep(l1outputfile,'tb','ta');
if exist(ncoutputfile); delete(ncoutputfile); end; % delete lv2 netcdf file if already exists

% Select variables to write
% Variables to be written in the lev2 netcdf file
vars = strvcat('time','time_bnds','lat','lon','zsl','azi','ele','height','ta','ta_offset','ta_err','ta_sys');

% Converting seconds from midnight to seconds from 1/1/1970
% NB: Be careful here: does O.time and R matches? How about non convergence episodes?
time_bnds = C.sampling;
time_offset = ( datenum(C.day_one(1),C.day_one(2),C.day_one(3)) - datenum(1970,1,1) ) * 24 * 3600;

% Interpolate on a fixed height grid
nret = length(R);
nlev = length(C.FixZgrid);
time = zeros(nret,1); 
ta     = C.MissingValue * ones(nlev,nret); % set to missing value
ta_err = C.MissingValue * ones(nlev,nret); % set to missing value
ta_sys = C.MissingValue * ones(nlev,nret); % set to missing value
for ir = 1:nret
    time(ir) = R(ir).time + time_offset;
    if R(ir).nite <= C.MaxIterations % convergence was reached
       ibkg = X.ObsMatch(R(ir).nobs_prf);       
       if ( X.Z(1,ibkg) < C.FixZgrid(1) ) | ( X.Z(end,ibkg) > C.FixZgrid(end) ); 
          disp('Interpolation may produce NaN at the vertical boudaries'); return; 
       end
       ta(:,ir)     = interp1( X.Z(:,ibkg) , R(ir).Ret_T_K , C.FixZgrid);
       ta_err(:,ir) = interp1( X.Z(:,ibkg) , A(ir).TTotErr , C.FixZgrid);
       ta_sys(:,ir) = interp1( X.Z(:,ibkg) , E(ir).TSysUnc , C.FixZgrid); 
    end
end
ta_offset = zeros(size(ta_err(:,1))); % I leave this field zero unless it's needed in the future

% Create variables  
nccreate(ncoutputfile,'time','Dimensions',{'time' length(time)},'Datatype','double','Format','netcdf4');
nccreate(ncoutputfile,'time_bnds','Dimensions',{'nv' 2 'time'},'Datatype','double');
nccreate(ncoutputfile,'lat','Datatype','single');
nccreate(ncoutputfile,'lon','Datatype','single');
nccreate(ncoutputfile,'zsl','Datatype','single');
nccreate(ncoutputfile,'azi','Dimensions',{'n_angle' length(O.angles_az)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'ele','Dimensions',{'n_angle' length(O.angles_el)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'height','Dimensions',{'height' nlev},'Datatype','single');
nccreate(ncoutputfile,'ta','Dimensions',{'height' 'time'},'Datatype','single');  
nccreate(ncoutputfile,'ta_offset','Dimensions',{'height'},'Datatype','single');      % different than original lv2 (ta_offset has height and time dimension)
nccreate(ncoutputfile,'ta_sys','Dimensions',{'height' 'time'},'Datatype','single');  % new wrt to original lv2. ta_sys has height and time dimension
nccreate(ncoutputfile,'ta_err','Dimensions',{'height' 'time'},'Datatype','single');  % different than original lv2 (ta_err had only height dimension)


% Write data associated to variables
ncwrite(ncoutputfile,'time',time);
ncwrite(ncoutputfile,'time_bnds',[time-time_bnds time+time_bnds]');
ncwrite(ncoutputfile,'lat',C.lat);
ncwrite(ncoutputfile,'lon',C.lon);
ncwrite(ncoutputfile,'zsl',C.asl);
ncwrite(ncoutputfile,'azi',O.angles_az);
ncwrite(ncoutputfile,'ele',O.angles_el);
ncwrite(ncoutputfile,'height',C.FixZgrid);
ncwrite(ncoutputfile,'ta',ta); 
ncwrite(ncoutputfile,'ta_offset',ta_offset); 
ncwrite(ncoutputfile,'ta_err',ta_err);
ncwrite(ncoutputfile,'ta_sys',ta_sys);

% Set variable attributes
Variables = Net1DVariables_lv2;
nvar = length(Variables);
for iv = 1:nvar
    varName = Variables(iv).Name;  
    if ismember(varName,vars,'rows')        
              ncwriteallatt(ncoutputfile,Variables(iv));
    end    
end
    
% Set global attributes
ncwriteatt(ncoutputfile,'/','Title','Microwave radiometer retrieved temperature profile');
ncwriteatt(ncoutputfile,'/','Institution',[O.GA.Institution]);
ncwriteatt(ncoutputfile,'/','Contact_person',[C.nc.glatt.Contact_person O.GA.Contact_person]);
ncwriteatt(ncoutputfile,'/','Source',[O.GA.Source]);
ncwriteatt(ncoutputfile,'/','History',[C.nc.glatt.History]);
ncwriteatt(ncoutputfile,'/','Dependencies',l1outputfile);
ncwriteatt(ncoutputfile,'/','Conventions',[O.GA.Conventions]);
ncwriteatt(ncoutputfile,'/','Processing_date',[C.nc.glatt.Processing_date]);
ncwriteatt(ncoutputfile,'/','Author',[C.nc.glatt.Author O.GA.Author]);
ncwriteatt(ncoutputfile,'/','Comments',[C.nc.glatt.Comments O.GA.Comments]);
ncwriteatt(ncoutputfile,'/','License',[O.GA.License]);
ncwriteatt(ncoutputfile,'/','Measurement_site',[O.GA.Measurement_site]);

return

