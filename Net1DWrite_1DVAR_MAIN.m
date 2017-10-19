%main function to write the input files to the NWPSAF 1DVAR

function X_new=Net1DWrite_1DVAR_MAIN(O,X,C)

% O = structure with observation information
% X = structure with the background information
% C = structure with configuration of the 1DVar 

% OPEN OUTPUT FILES FOR 1DVAR
fid_back=fopen([C.ODVARpath,'Sample_MWR/BACKGROUND.dat'],'w');
fid_obs=fopen([C.ODVARpath,'Sample_MWR/ObsFile.dat'],'w');

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
        %Qsurf=qspec_to_ppmv(X.Q(end,idx_back));  % Modified after checking with Pauline (see email of 2017/10/19) 
        Qspec=qabsolue_to_qspec(X.P(:,idx_back),X.T(:,idx_back),X.Q(:,idx_back));
        Qsurf=qspec_to_ppmv(Qspec(end));          % Modified after checking with Pauline (see email of 2017/10/19)
        
        % definition of the input CLW profile depending on the
        % configuration chosen
         % default value for C.Skyconditions==1, background put to zero and large background errors
        
        if O.cld31(t)==1
            CLW_tempo=zeros(X.nlev,1);
            if C.Skyconditions==2  % background LWP scaled by observed LWP
                if X.LWP(idx_back)~=0
                    CLW_tempo=X.LWC(:,idx_back)*O.LWP(t)/X.LWP(idx_back);
                else                           
                    r=Qspec./(1-Qspec);       % rapport de m??lange en kg/kg
                    es=6.107*10.^((7.5*(X.T(:,idx_back)-273.15))./(237.3+(X.T(:,idx_back)-273.15))); % Formule de T??tens donnant la pression de vapeur saturante en hPa
                    e=(X.P(:,idx_back))./(1+(0.622./r)); % pression partielle de vapeur d??duite de la pression et du rapport de m??lange en Pa
                    HR=100*e./(es*100); % Humidit?? relative en %                                       
               %    level_maxQ=find(X.Q(:,idx_back)==max(X.Q(:,idx_back)));
                   CLW_tempo(HR==max(HR))=1e-8;
                   LWP_tempo=0;
                    for lev=1:X.nlev-1
                        LWP_tempo=LWP_tempo+(CLW_tempo(lev)+CLW_tempo(lev+1))/2 *(X.P(lev+1,idx_back)-X.P(lev,idx_back))/9.81;
                    end          
 
                   CLW_tempo=CLW_tempo*O.LWP(t)/LWP_tempo;
                end
            elseif C.Skyconditions==3; % we follow the background profile
               CLW_tempo=X.LWC(:,idx_back);
              if X.LWP(idx_back)<1e-8
                   r=Qspec./(1-Qspec);       % rapport de m??lange en kg/kg
                   es=6.107*10.^((7.5*(X.T(:,idx_back)-273.15))./(237.3+(X.T(:,idx_back)-273.15))); % Formule de T??tens donnant la pression de vapeur saturante en hPa
                   e=(X.P(:,idx_back))./(1+(0.622./r)); % pression partielle de vapeur d??duite de la pression et du rapport de m??lange en Pa
                   HR=100*e./(es*100); % Humidit?? relative en %   
                   CLW_tempo(HR==max(HR))=1e-8;
              end
            end
        else
            CLW_tempo=X.LWC(:,idx_back);
        end
               
        % input Background file for 1DVAR
        Net1DWrite_BACK_1DVAR(X.P(:,idx_back),X.T(:,idx_back),Qspec(:),CLW_tempo,X.ps(idx_back),Qsurf,fid_back,nb_obs);

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
