function Net1DWrite_header_BACK_1DVAR(nlev,nb_prof,fid)

% function to write the header of the background file used for the 1DVAR retrieval
% inputs: C structure from the configuration file
% nlev : number of vertical levels of the background profiles
% nb_prof: total number of profiles that will be written in the background file
% fid = file identifier of the output file

    fprintf(fid,'%s \n','Background file generated from  ');
    fprintf(fid,'%s \n','AROME MODEL');
    fprintf(fid,'%s \n',sprintf('%d fixed layers',nlev));
    fprintf(fid,'%s \n','Created by P.Martinet');
    fprintf(fid,'%s \n','                           ');
    fprintf(fid,'%s \n','                           ');
    fprintf(fid,'%s \n','                           ');
    fprintf(fid,'%s \n','                           ');
    fprintf(fid,'%s \n','                           ');
    fprintf(fid,'%s \n','x------------------- End of Header------------------------------------x ');
    fprintf(fid,'%s %d \n','No. Background Profiles:  ',nb_prof);
    fprintf(fid,'%s \n',sprintf('No. of Levels/Profile:             %d',nlev));
    fprintf(fid,'%s \n','Unit for q:                         2  (1=ppmv, 2=kg/kg, 3=RH) ');
