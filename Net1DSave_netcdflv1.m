% Network 1DVAR+RTTOV retrieval: Save 1DVAR output data in netcdf
%
% Net1DSave_netcdflv1 saves output level1 data according to Config C structure
% The output file is a CF-compliant netcdf file, as checked through: 
% http://puma.nerc.ac.uk/cgi-bin/cf-checker.pl

function Net1DSave_netcdflv1(ncoutputfile,C,O);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select variables to write
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables to be written in the lev1 netcdf file
vars = strvcat('time','time_bnds','lat','lon','zsl','freq_sb','azi','ele','tb','offset_tb','freq_shift','tb_bias','tb_cov','ta','pa','hur');

% Converting seconds from midnight to seconds from 1/1/1970
time_bnds = 60;
time_offset = ( datenum(C.day_one(1),C.day_one(2),C.day_one(3)) - datenum(1970,1,1) ) * 24 * 3600;
time = O.time + time_offset; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create variables  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nccreate(ncoutputfile,'time','Dimensions',{'time' length(O.time)},'Datatype','double','Format','netcdf4');
nccreate(ncoutputfile,'time_bnds','Dimensions',{'nv' 2 'time'},'Datatype','double');
nccreate(ncoutputfile,'lat','Datatype','single');
nccreate(ncoutputfile,'lon','Datatype','single');
nccreate(ncoutputfile,'zsl','Datatype','single');
nccreate(ncoutputfile,'freq_sb','Dimensions',{'n_freq' length(O.channels)},'Datatype','single');
nccreate(ncoutputfile,'azi','Dimensions',{'n_angle' length(O.angles_az)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'ele','Dimensions',{'n_angle' length(O.angles_el)},'Datatype','single');    % this may have to change for RAO!
nccreate(ncoutputfile,'tb','Dimensions',{'n_freq' 'n_angle' 'time'},'Datatype','single');         % this may have to change for RAO!
nccreate(ncoutputfile,'offset_tb','Dimensions',{'n_freq' 'n_angle' 'time'},'Datatype','single');  % this may have to change for RAO!
nccreate(ncoutputfile,'freq_shift','Dimensions',{'n_freq'},'Datatype','single');
nccreate(ncoutputfile,'tb_bias','Dimensions',{'n_freq'},'Datatype','single');
%nccreate(ncoutputfile,'tb_cov','Dimensions',{'n_freq' 'n_freq'},'Datatype','single');
%nccreate(ncoutputfile,'tb_cov','Dimensions',{'n_freq' 'n_freq' length(O.channels)},'Datatype','single');
nccreate(ncoutputfile,'tb_cov','Dimensions',{'n_freq' 'n_freq2' length(O.channels)},'Datatype','single');
nccreate(ncoutputfile,'ta','Dimensions',{'time'},'Datatype','single');
nccreate(ncoutputfile,'pa','Dimensions',{'time'},'Datatype','single');
nccreate(ncoutputfile,'hur','Dimensions',{'time'},'Datatype','single');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write data associated to variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncwrite(ncoutputfile,'time',time);
ncwrite(ncoutputfile,'time_bnds',[time-time_bnds time+time_bnds]');
ncwrite(ncoutputfile,'lat',C.lat);
ncwrite(ncoutputfile,'lon',C.lon);
ncwrite(ncoutputfile,'zsl',C.asl);
ncwrite(ncoutputfile,'freq_sb',O.channels);
ncwrite(ncoutputfile,'azi',O.angles_az);
ncwrite(ncoutputfile,'ele',O.angles_el);
ncwrite(ncoutputfile,'tb',O.y);
ncwrite(ncoutputfile,'offset_tb',O.ybias);
ncwrite(ncoutputfile,'freq_shift',O.freq_shift);
ncwrite(ncoutputfile,'tb_bias',sqrt(diag(O.Rcov))); % FixMe! Is this what should be?
ncwrite(ncoutputfile,'tb_cov',O.Rcov);
ncwrite(ncoutputfile,'ta',O.ta);
ncwrite(ncoutputfile,'pa',O.pa);
ncwrite(ncoutputfile,'hur',O.hur);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set variable attributes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nvar = length(O.ncinfo.Variables);
for iv = 1:nvar
    varName = O.ncinfo.Variables(iv).Name;  
    if ismember(varName,vars,'rows')        
              ncwriteallatt(ncoutputfile,O.ncinfo.Variables(iv));
    end    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set global attributes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncwriteatt(ncoutputfile,'/','Title',[O.GA.Title]);
ncwriteatt(ncoutputfile,'/','Institution',[O.GA.Institution]);
ncwriteatt(ncoutputfile,'/','Contact_person',[C.nc.glatt.Contact_person O.GA.Contact_person]);
ncwriteatt(ncoutputfile,'/','Source',[O.GA.Source]);
ncwriteatt(ncoutputfile,'/','History',[C.nc.glatt.History O.GA.History]);
ncwriteatt(ncoutputfile,'/','Conventions',[O.GA.Conventions]);
ncwriteatt(ncoutputfile,'/','Processing_date',[C.nc.glatt.Processing_date]);
ncwriteatt(ncoutputfile,'/','Author',[C.nc.glatt.Author O.GA.Author]);
ncwriteatt(ncoutputfile,'/','Comments',[C.nc.glatt.Comments O.GA.Comments]);
ncwriteatt(ncoutputfile,'/','License',[O.GA.License]);
ncwriteatt(ncoutputfile,'/','Measurement_site',[O.GA.Measurement_site]);
                              
return
