function Net1DWrite_OBS_1DVAR(TB,C,fid,obs_nb)

% function to write input observation file to the NWPSAF 1DVAR
% TB = matrix of observations; (size = number of instrument channel , number of elevation angle)
% C = Configuration structure from the config file
% fid= file identifier to write the data
% obs_nb = number/position of the observation to be process among the whole dataset of observations

    fprintf(fid,'%s \n',sprintf('Obs ID:              %d Obs Type:          3 Satellite ID:   1 ',obs_nb));
    fprintf(fid,'%s \n',sprintf(' Latitude:  %.2f Longitude:   %2f Elevation:    %f ',C.lat,C.lon,C.asl));
    fprintf(fid,'%s \n',' Surface Type:   1 Sat Zen Angle:    0 Solar Zen. Ang.:    0.000  ');    
    fprintf(fid,'%s \n','Brightness temperatures');
    for ang=1:C.elev_angles_how_many
        fprintf(fid,'%f \n',TB(C.elev_angles_channum{ang},ang));
    end
       
end
