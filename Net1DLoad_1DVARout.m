% Network 1DVAR+RTTOV retrieval: Load 1DVAR output data
%
% Net1DLoad_1DVARout loads output 1DVAR data according to Config C 
% structure, i.e. for different station, day,...

function R = Net1DLoad_1DVARout(C,X);

% initialize R structure
empty = [];
R = struct('nobs_prf',empty,'nobs_tbs',empty,'Bkg_P_hPa',empty,'Bkg_T_K',empty,'Bkg_Q_kgkg',empty,'Ret_T_K',empty,'Ret_Q_kgkg',empty,'nchn',empty,'Bkg_Tb_K',empty,'Ret_Tb_K',empty,'Obs_Tb_K',empty);


% find output files
out_path = C.ODVARpath;
out_prof = [out_path '/Retrieved_Profiles.dat'];
out_tbs =  [out_path '/Retrieved_BTs.dat'];


% open output files
fid_prof = fopen(out_prof,'r');
fid_tbs  = fopen(out_tbs,'r');


% Loop on output profile and tbs file
npro = 0;
while 1
    
    line = fgetl(fid_prof);
    if strcmp(line(1:12),' Observation'); % here start a new retrieved profile

        npro = npro + 1;
        nobs = str2num(line(20:27));

        % skip header lines
        line = fgetl(fid_prof);
        line = fgetl(fid_prof);
        MTX = fscanf(fid_prof,'%g',[7,X.nlev]); MTX = MTX';
        line = fgetl(fid_prof);
        line = fgetl(fid_prof);
        line = fgetl(fid_prof);
        if C.retrieve_LWP==1
            line = fgetl(fid_prof);
            nite = str2num(line(30:end));
        else
            nite = str2num(line(30:end));
        end
        
        % fill the R structure
        R(npro).nite       = nite;
        R(npro).nobs_prf   = nobs;
        R(npro).Bkg_P_hPa  = MTX(:,1);
        R(npro).Bkg_T_K    = MTX(:,5);
        R(npro).Bkg_Q_kgkg = MTX(:,6);
        if C.retrieve_T(1)
        R(npro).Ret_T_K    = MTX(:,2);
        end
        if C.retrieve_Q(1)
        R(npro).Ret_Q_kgkg = MTX(:,3);
        end
       
    end
    
    if feof(fid_prof); break; end;
    
end

% Loop on output Tb file
npro = 0;
while 1
    
    line = fgetl(fid_tbs);
    if strcmp(line(1:12),' Observation'); % here start a new obs

        npro = npro + 1;
        nobs = str2num(line(20:27));
        line = fgetl(fid_tbs);
        nchn = str2num(line(30:39));

        % skip 1 header line1
        line = fgetl(fid_tbs);
        MTX = fscanf(fid_tbs,'%g',[4,nchn]); MTX = MTX';
        line = fgetl(fid_tbs);
        
        % fill the R structure
        R(npro).nobs_tbs = nobs;
        R(npro).nchn     = nchn;
        R(npro).Bkg_Tb_K = MTX(:,2);
        R(npro).Ret_Tb_K = MTX(:,4);
        R(npro).Obs_Tb_K = MTX(:,3);
        
    end
    
    if feof(fid_tbs); break; end;
    
end


% close output files
fclose(fid_prof);
fclose(fid_tbs);

return


