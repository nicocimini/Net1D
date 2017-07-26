% Network 1DVAR+RTTOV retrieval: Save 1DVAR output data
%
% Net1DLoad_save1DVARout saves output 1DVAR data according to Config C 
% structure, i.e. for different station, day,...

function Net1DSave_1DVARout(C,O,X,R,A,E);

% Pauline: save original output 1DVAR files
out_path = [C.ODVARpath_output C.station_id '/'];

if ~exist(out_path,'dir')
    mkdir(out_path)
end

if C.GeneralMode_1DVAR < 30
    liste_files={'Retrieved_BTs.dat','Retrieved_Profiles.dat'};
else
    liste_files={'Retrieved_BTs.dat','Retrieved_Profiles.dat','Minimisation_BT.log',...
        'Minimisation.log','ProfileQC.dat','A-Matrix.out','Am-Matrix.out',...
        'AveragingKernel.out','BgJacobian.out','RetJacobian.out'};
end

for i=1:length(liste_files)
   fname_output=[liste_files{i}(1:end-4) '_' C.station_id '_' num2str(C.day_one(1)) twodigstr(C.day_one(2)) twodigstr(C.day_one(3))  liste_files{i}(end-3:end)];
   command=['mv ' C.ODVARpath_retrieval liste_files{i} ' ' out_path fname_output  ];
   [status,result] = system(command);
end


switch C.Output_format 
      
    case 'matlab'; 
        
          % naming output path and files
          out_path = [C.ODVARpath_output C.station_id '/'];
          out_file = ['Net1D_output_' num2str(C.day_one(1)) twodigstr(C.day_one(2)) twodigstr(C.day_one(3))];
          save([out_path out_file '.mat'],'C','X','R','A','E');
        
    case 'netcdf'

          % lev1
          indx = strfind(O.ncinfo.Filename,'/');
          out_path_lv1 = [C.NetCDF_output 'level1/' C.station_id '/'];
          out_file_lv1 = O.ncinfo.Filename(indx(end)+1:end);
          %Net1DSave_netcdflv1([out_path_lv1 out_file_lv1],C,O);
          
          % lev2
          out_path_lv2 = strrep(out_path_lv1,'level1','level2');
          out_file_lv2 = strrep(out_file_lv1,'l1','l2');
          % root filename changes inside for different output files (ta, hua, ...)
          Net1DSave_netcdflv2([out_path_lv2 out_file_lv2],C,O,X,R,A);
        
end

return