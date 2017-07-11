%main function to write the input files to the NWPSAF 1DVAR

function X_new=Net1DWrite_1DVAR_MAIN(O,X,C)

% O = structure with observation information
% X = structure with the background information
% C = structure with configuration of the 1DVar 

% OPEN OUTPUT FILES FOR 1DVAR
fid_back=fopen([C.ODVARpath,'BACKGROUND.dat'],'w');
fid_obs=fopen([C.ODVARpath,'ObsFile.dat'],'w');

% how many observations to be process
nb_obs=0;

%Index of the background used for each profile in the 1DVAR
X.ObsMatch=zeros(length(O.time),1);
X.ObsMatch=-9999; %default value if the observation has not been selected
% Flag 1/0 observation selected for the 1DVAR (clear-sky/background timelapse selected)
X.ObsUsed=zeros(length(O.time),1);

for t=1:length(O.time) % loop over all the observations
%for t=1:1
    % find background profile closest to observation
    idx_back=find(X.time<O.time(t),1,'last'); 
    deltatime_obs_background=abs(X.time(idx_back)-O.time(t));
    % continue only if background less than 6 hours from observation
    if deltatime_obs_background < 6*3600
        nb_obs=nb_obs+1;
    end
end

%write header of the observation file 
Net1DWrite_header_OBS_1DVAR(C,nb_obs,fid_obs)

%write header of the background file 
Net1DWrite_header_BACK_1DVAR(X.nlev,nb_obs,fid_back)

nb_obs=0;
for t=1:length(O.time) % loop over all the observations
%for t=1:1
    
    % find background profile closest to observation
    idx_back=find(X.time<O.time(t),1,'last'); 
    deltatime_obs_background=abs(X.time(idx_back)-O.time(t));
    % continue only if background less than 6 hours from observation
    
    if deltatime_obs_background < 6*3600
        X.ObsUsed(t)=1; % observation can be used for 1DVAR
        X.ObsMatch(t)=idx_back;
        nb_obs=nb_obs+1;
        %computation of surface humidity in ppmv for 1DVAR
        Qsurf=qspec_to_ppmv(X.Q(end,idx_back));
        Qspec=qabsolue_to_qspec(X.P(:,idx_back),X.T(:,idx_back),X.Q(:,idx_back));
        
        % input Background file for 1DVAR
        Net1DWrite_BACK_1DVAR(X.P(:,idx_back),X.T(:,idx_back),Qspec(:),X.LWC(:,idx_back),X.ps(idx_back),Qsurf,fid_back,nb_obs);

        % input Obs file for 1DVAR
        Net1DWrite_OBS_1DVAR(squeeze(O.y(:,:,t)),C,fid_obs,nb_obs);
    end        
end

%input R matrix for 1DVAR
Net1DWrite_Rmatrix_1DVAR(C,O.Rcov)

%input B matrix for 1DVAR
Net1DWrite_Bmatrix_1DVAR(X,C)

%input 1DVAR control file
Net1DWrite_Controlfile_1DVAR(C)

%input 1DVAR retrieval file
Net1DWrite_Retrievalfile_1DVAR(C, X.nlev)

% input 1DVAR channelchoice
Net1DWrite_Channelchoice_1DVAR(C)

X_new=X;

fclose(fid_back);
fclose(fid_obs);
