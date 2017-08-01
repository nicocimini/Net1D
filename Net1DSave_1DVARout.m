% Network 1DVAR+RTTOV retrieval: Save 1DVAR output data
%
% Net1DLoad_save1DVARout saves output 1DVAR data according to Config C 
% structure, i.e. for different station, day,...

function Net1DSave_1DVARout(C,O,X,R,A,E);

% Save output files according to C.Output_format_numb


% Loop on output formats
for numb = C.Output_format_numb
    
    Output_format = C.Output_format_name(numb,:);
    % Create folder if not existing
    out_path = [C.ODVARpath_output Output_format '/' C.station_id '/'];
    if ~exist(out_path,'dir')
       mkdir(out_path)
    end

    switch Output_format 
 
    case 'NWPSAF'; % save original output 1DVAR files
        
        % Move ascii files from original dir to output dir
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
              
    case 'matlab'; % save in matlab format
        
          % naming output path and files
          %out_path = [C.ODVARpath_output C.station_id '/'];
          out_file = ['Net1D_output_' num2str(C.day_one(1)) twodigstr(C.day_one(2)) twodigstr(C.day_one(3))];
          save([out_path out_file '.mat'],'C','X','R','A','E');
        
    case 'netcdf' % save in HDCP2 netcdf format

          % lev1
          indx = strfind(O.ncinfo.Filename,'/');
          %out_path_lv1 = [C.NetCDF_output 'level1/' C.station_id '/'];
          out_path_lv1 = [out_path 'level1/']; 
          if ~exist(out_path_lv1,'dir'); mkdir(out_path_lv1); end
          out_file_lv1 = O.ncinfo.Filename(indx(end)+1:end);
          Net1DSave_netcdflv1([out_path_lv1 out_file_lv1],C,O);
          
          % lev2
          out_path_lv2 = strrep(out_path_lv1,'level1','level2');
          if ~exist(out_path_lv2,'dir'); mkdir(out_path_lv2); end
          out_file_lv2 = strrep(out_file_lv1,'l1','l2');
          % root filename changes inside for different output files (ta, hua, ...)
          Net1DSave_netcdflv2([out_path_lv2 out_file_lv2],C,O,X,R,A);
        
    end
    
end

return