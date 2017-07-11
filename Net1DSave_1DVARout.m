% Network 1DVAR+RTTOV retrieval: Save 1DVAR output data
%
% Net1DLoad_save1DVARout saves output 1DVAR data according to Config C 
% structure, i.e. for different station, day,...

function Net1DSave_1DVARout(C,O,X,R,A,E);

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