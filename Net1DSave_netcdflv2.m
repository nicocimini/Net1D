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

function Net1DSave_netcdflv2(l1outputfile,C,O,X,R,A);


% Write ta file
% lv2file = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/level2/joy/sups_joy_mwrBL00_l2_ta_v01_20140101005346.nc'
% ncdisp(lv2file);
if C.retrieve_T(1)
   Net1DSave_netcdflv2_ta(l1outputfile,C,O,X,R,A);
end

% Write hua file
% lv2file = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/level2/joy/sups_joy_mwr00_l2_hua_i00_20140101000004.nc'
% ncdisp(lv2file);
if C.retrieve_Q(1)
   Net1DSave_netcdflv2_hua(l1outputfile,C,O,X,R,A);
end

% Write prw file
% lv2file = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/level2/joy/sups_joy_mwr00_l2_prw_i00_20140101000004.nc'
% ncdisp(lv2file);

% Write clwvi file
% lv2file = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/level2/joy/sups_joy_mwr00_l2_clwvi_i00_20140101000004.nc'
% ncdisp(lv2file);


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub-functions to write different NetCDF file types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Net1DSave_netcdflv2_hua %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Net1DSave_netcdflv2_hua(l1outputfile,C,O,X,R,A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change output name (to do for each variable/file, e.g. hua, ipwv,...)
ncoutputfile = strrep(l1outputfile,'tb','hua');

% Select variables to write
% Variables to be written in the lev2 netcdf file
vars = strvcat('time','time_bnds','lat','lon','zsl','azi','ele','ele_ret','height','hua','hua_offset','hua_err');

% Converting seconds from midnight to seconds from 1/1/1970
% NB: Be careful here: does O.time and R matches? How about non convergence episodes? 
% If R(i).nite >= C.MaxIterations means no convergence. Should treat this!
time_bnds = 60;
time_offset = ( datenum(C.day_one(1),C.day_one(2),C.day_one(3)) - datenum(1970,1,1) ) * 24 * 3600;
time = O.time + time_offset; 

% Interpolate on a fixed height grid
nret = length(R);
nlev = length(C.FixZgrid);
hua     = C.MissingValue * ones(nlev,nret); % set to missing value
hua_err = C.MissingValue * ones(nlev,nret); % set to missing value
for ir = 1:nret
    if R(ir).nite <= C.MaxIterations % convergence was reached
       % from specific (kg/kg) to absolute (kg/m-3) humidity
       Qa = qspec_to_qabsolue(R(ir).Bkg_P_hPa*100,R(ir).Ret_T_K,R(ir).Ret_Q_kgkg);
       Qa_err = qspec_to_qabsolue(R(ir).Bkg_P_hPa*100,R(ir).Ret_T_K,sqrt(diag(A(ir).AQ))); % NB: A of Q
       ibkg = X.ObsMatch(ir);
       if ( X.Z(1,ibkg) < C.FixZgrid(1) ) | ( X.Z(end,ibkg) > C.FixZgrid(end) ); 
          disp('Interpolation may produce NaN at the vertical boudaries'); return; 
       end
       hua(:,ir)     = interp1( X.Z(:,ibkg) , Qa     , C.FixZgrid);
       hua_err(:,ir) = interp1( X.Z(:,ibkg) , Qa_err , C.FixZgrid); 
    end
end
hua_offset = zeros(size(hua_err(:,1))); % I leave this field zero unless it's needed in the future


% Create variables  
nccreate(ncoutputfile,'time','Dimensions',{'time' length(O.time)},'Datatype','double','Format','netcdf4');
nccreate(ncoutputfile,'time_bnds','Dimensions',{'nv' 2 'time'},'Datatype','double');
nccreate(ncoutputfile,'lat','Datatype','single');
nccreate(ncoutputfile,'lon','Datatype','single');
nccreate(ncoutputfile,'zsl','Datatype','single');
nccreate(ncoutputfile,'azi','Dimensions',{'n_angle' length(O.angles_az)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'ele','Dimensions',{'n_angle' length(O.angles_el)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'ele_ret','Dimensions',{'time'},'Datatype','single');                       % this may have to change for RAO!
nccreate(ncoutputfile,'height','Dimensions',{'height' nlev},'Datatype','single');
nccreate(ncoutputfile,'hua','Dimensions',{'height' 'time'},'Datatype','single');  
nccreate(ncoutputfile,'hua_offset','Dimensions',{'height'},'Datatype','single');      % different than original lv2 (hua_offset has height and time dimension)
nccreate(ncoutputfile,'hua_err','Dimensions',{'height' 'time'},'Datatype','single');  % different than original lv2 (hua_err has only height dimension)

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
ncwrite(ncoutputfile,'hua_offset',hua_offset);  % we probably do not want to output an offset for ta... to be decided!
ncwrite(ncoutputfile,'hua_err',hua_err);

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
function Net1DSave_netcdflv2_ta(l1outputfile,C,O,X,R,A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Change output name (to do for each variable/file, e.g. hua, ipwv,...)
ncoutputfile = strrep(l1outputfile,'tb','ta');

% Select variables to write
% Variables to be written in the lev2 netcdf file
vars = strvcat('time','time_bnds','lat','lon','zsl','azi','ele','height','ta','ta_offset','ta_err');

% Converting seconds from midnight to seconds from 1/1/1970
% NB: Be careful here: does O.time and R matches? How about non convergence episodes?
time_bnds = 60;
time_offset = ( datenum(C.day_one(1),C.day_one(2),C.day_one(3)) - datenum(1970,1,1) ) * 24 * 3600;
time = O.time + time_offset; 

% Interpolate on a fixed height grid
nret = length(R);
nlev = length(C.FixZgrid);
ta     = C.MissingValue * ones(nlev,nret); % set to missing value
ta_err = C.MissingValue * ones(nlev,nret); % set to missing value
for ir = 1:nret
    if R(ir).nite <= C.MaxIterations % convergence was reached
       ibkg = X.ObsMatch(ir);
       if ( X.Z(1,ibkg) < C.FixZgrid(1) ) | ( X.Z(end,ibkg) > C.FixZgrid(end) ); 
          disp('Interpolation may produce NaN at the vertical boudaries'); return; 
       end
       ta(:,ir)     = interp1( X.Z(:,ibkg) , R(ir).Ret_T_K       , C.FixZgrid);
       ta_err(:,ir) = interp1( X.Z(:,ibkg) , sqrt(diag(A(ir).AT)) , C.FixZgrid);
    end
end
ta_offset = zeros(size(ta_err(:,1))); % I leave this field zero unless it's needed in the future

% Create variables  
nccreate(ncoutputfile,'time','Dimensions',{'time' length(O.time)},'Datatype','double','Format','netcdf4');
nccreate(ncoutputfile,'time_bnds','Dimensions',{'nv' 2 'time'},'Datatype','double');
nccreate(ncoutputfile,'lat','Datatype','single');
nccreate(ncoutputfile,'lon','Datatype','single');
nccreate(ncoutputfile,'zsl','Datatype','single');
nccreate(ncoutputfile,'azi','Dimensions',{'n_angle' length(O.angles_az)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'ele','Dimensions',{'n_angle' length(O.angles_el)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'height','Dimensions',{'height' nlev},'Datatype','single');
nccreate(ncoutputfile,'ta','Dimensions',{'height' 'time'},'Datatype','single');  
nccreate(ncoutputfile,'ta_offset','Dimensions',{'height'},'Datatype','single');      % different than original lv2 (ta_offset has height and time dimension)
nccreate(ncoutputfile,'ta_err','Dimensions',{'height' 'time'},'Datatype','single');  % different than original lv2 (ta_err has only height dimension)


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
ncwrite(ncoutputfile,'ta_offset',ta_offset);  % we probably do not want to output an offset for ta... to be decided!
ncwrite(ncoutputfile,'ta_err',ta_err);

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

