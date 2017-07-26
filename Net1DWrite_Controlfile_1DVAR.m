function Net1DWrite_Controlfile_1DVAR(C)

% function to write the 1DVAR input control file
% input: C structure with 1DVAR init configuration

fid=fopen([C.ODVARpath,'Sample_MWR/','ControlData.NL'],'w');

if strcmp(C.instrument,'HATPRO')
    fprintf(fid,'%s \n','! HATPRO ControlFile for NWPSAF 1DVAR '); 
    COEFFS_DIR='HATPRO_COEFFS_DIR/';
elseif strcmp(C.instrument,'MP3000A')
    fprintf(fid,'%s \n','! MP3000A ControlFile for NWPSAF 1DVAR '); 
    COEFFS_DIR='MP3000A_COEFFS_DIR/';
else
    disp('WRONG INSTRUMENT')
end

if strcmp(C.Minimisation_Method,'ML')
    Minimisation_Method=2;
elseif strcmp(C.Minimisation_Method,'N')
    Minimisation_Method=1;
else
    disp('WRONG MINIMISATION METHOD')
end

fprintf(fid,'%s \n',''); 
fprintf(fid,'%s \n','! ******** THESE COMMENTS MAY NEED TO BE REMOVED IF COMPILING ***********'); 
fprintf(fid,'%s \n','! ******** WITH F90 (rather than F95) ***********************************'); 
fprintf(fid,'%s \n','');
%fprintf(fid,'%s \n','! Relative path of the coefficients directory');
fprintf(fid,'%s \n','&Control');
fprintf(fid,'%s \n','! This controls the verbosity of the output 0=Minimal, 40=verbose '); 
fprintf(fid,'%s %d \n','GeneralMode=',C.GeneralMode_1DVAR);  
fprintf(fid,'%s \n','! For the minimisation option, 1 is Newtonian, 2 is Marquardt-Levenberg'); 
fprintf(fid,'%s %d \n','Minimisation_Method=',Minimisation_Method); 
fprintf(fid,'%s %d \n','MaxIterations=',C.MaxIterations); 
fprintf(fid,'%s \n','Legacy_Setting=.true.'); 
%fprintf(fid,'%s \n','! Cloud Liquid Water profile   '); 
fprintf(fid,'%s \n','! These are cloud cost thresholds:'); 
fprintf(fid,'%s \n','CostThresh_Land=20. '); 
fprintf(fid,'%s \n','CostThresh_Sea=20.'); 
fprintf(fid,'%s \n','CostThresh_IRWindow_Land=0.'); 
fprintf(fid,'%s \n','CostThresh_IRWindow_Sea=0.');
fprintf(fid,'%s \n','CloudAbsThresh_IRWindow=-10. '); 
fprintf(fid,'%s \n','HighCloudAbsThresh_IRWindow=-55.');
fprintf(fid,'%s \n','retrieval_in_log=.false.');
if C.force_bgk_clear ==1
    fprintf(fid,'%s \n','Read_CLW_Background =.false.');
else
    fprintf(fid,'%s \n','Read_CLW_Background =.true.');
end

fprintf(fid,'%s %f %s \n','LwpSD=',C.LWPsd,'/');
fprintf(fid,'%s \n','       ');
fclose(fid);






