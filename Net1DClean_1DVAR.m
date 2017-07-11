%function to clean the 1DVAR directory of previous output files
% + link to the right RTTOV directory

function Net1DClean_1DVAR(C)

% find output files
out_path =  C.ODVARpath;
out_prof = [out_path '/Retrieved_Profiles.dat'];
out_tbs =  [out_path '/Retrieved_BTs.dat'];
out_AK= [out_path '/AveragingKernel.out'];
out_Jac=[out_path '/RetJacobian.out'];
out_A = [out_path '/A-Matrix.out'];
out_Am = [out_path '/Am-Matrix.out'];
out_mini= [out_path '/Minimisation.log'];
out_mini_TB= [out_path '/Minimisation_BT.log'];

%input files to be removed also
fid_back=[C.ODVARpath,'BACKGROUND.dat'];
fid_obs=[C.ODVARpath,'ObsFile.dat'];
fid_retrieval=[C.ODVARpath,'Retrieval.NL'];
fid_controlfile= [C.ODVARpath,'ControlData.NL'];

if strcmp(C.instrument,'HATPRO')
    fid_B=[C.ODVARpath,'HATPRO_COEFFS_DIR/','Bmatrix'];
    fid_R= [C.ODVARpath,'HATPRO_COEFFS_DIR/','Rmatrix'];
    fid_channel_choice= [C.ODVARpath,'HATPRO_COEFFS_DIR/','ChannelChoice.dat'];
elseif strcmp(C.instrument,'MP3000A')  
    fid_B=[C.ODVARpath,'MP3000A_COEFFS_DIR/','Bmatrix'];
    fid_R= [C.ODVARpath,'MP3000A_COEFFS_DIR/','Rmatrix'];
    fid_channel_choice= [C.ODVARpath,'MP3000A_COEFFS_DIR/','ChannelChoice.dat'];        
end

command=['rm ' out_prof ' ' out_prof ' ' out_tbs ' ' out_AK ' ' out_Jac ' ' out_A ' ' out_Am ' ' out_mini ' ' out_mini_TB];
[status,result] = system(command);
command=['rm ' fid_back ' ' fid_obs ' ' fid_retrieval ' ' fid_controlfile];
[status,result] = system(command);
command=['rm ' fid_B ' ' fid_R ' ' fid_channel_choice];
[status,result] = system(command);

%symbolink link to the good rttov coefficient
