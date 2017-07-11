function Net1DWrite_header_OBS_1DVAR(C,nb_obs,fid)

% function to write the header of the observation file used for the 1DVAR retrieval
% inputs: C structure from the configuration file
% nb_obs: total number of observations that will be written in the observation file
% fid = file identifier of the output file

nb_chan=0;
for ang=1:C.elev_angles_how_many
    nb_chan=nb_chan+length(C.elev_angles_channum{ang});
end

if strcmp(C.instrument,'HATPRO')
    instrument_number=75;
    fprintf(fid,'%s \n',' HATPRO OBSERVATION '); 
elseif strcmp(C.instrument,'MP3000A')
    instrument_number=76;
    fprintf(fid,'%s \n',' MP3000A OBSERVATION '); 
end
  
    fprintf(fid,'%s %s\n','Created',date);
    for i=1:8
        fprintf(fid,'%s \n','                           ');
    end
    fprintf(fid,'%s %d \n','Number of Observations in File:  ', nb_obs); 
    fprintf(fid,'%s \n',sprintf('No. of Chans per Observation:      %d',nb_chan));
    fprintf(fid,'%s \n',sprintf('Number of instruments making up observations : %d ',C.elev_angles_how_many));
    fprintf(fid,'%s \n','*** In the following Series, Platform and Instrument are defined  *** ');
    fprintf(fid,'%s \n','*** according to the relevant RT Model definitions (if required): *** ');
    fprintf(fid,'%s \n','Composite Instruments: 1');  
    if strcmp(C.instrument,'HATPRO')
        fprintf(fid,'%s \n','HATPRO'); 
    elseif strcmp(C.instrument,'MP3000A')
        fprintf(fid,'%s \n','MP3000A'); 
    end
    fprintf(fid,'%s \n','Sat. Series   Platform   Instrument First_Channel   Last_Channel  Sat ID View Angle');
    first_chan=1;
    for ang=1:C.elev_angles_how_many
        last_chan=first_chan+length(C.elev_angles_channum{ang})-1;
        fprintf(fid,'%s \n',sprintf('   41        1      %d       %d       %d     1     %4.2f ',instrument_number,first_chan,last_chan,90-C.elev_angles_degrees(ang)));
        first_chan=last_chan+1;
    end
    fprintf(fid,'%s \n','Channels:'); 
    for ang=1:C.elev_angles_how_many
        fprintf(fid,'%s',sprintf('  %d',C.elev_angles_channum{ang})); 
    end
    fprintf(fid,'%s \n','    '); %%SAUT DE LIGNE
    fprintf(fid,'%s \n','------------------------------------------------------- ');
