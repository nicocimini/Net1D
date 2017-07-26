% Net1DLoad_1DVARout loads output 1DVAR data according to Config C 
% structure, i.e. for different station, day,...
% only averaging kernal matrix

function R = Net1DLoad_1DVARout_AK(C,X)

% find output files

%out_path = [C.ODVARpath C.station_id];
out_AK = [C.ODVARpath_retrieval 'AveragingKernel.out'];


% open output files
fid_AK = fopen(out_AK,'r');

nb_prof=0;
while 1
    
    line = fgetl(fid_AK);
  % disp(line)
    if strcmp(line(1:9),' ObNumber'); % here start a new retrieved profile
     %   nobs = str2num(line(22:27));
        line = fgetl(fid_AK);
        %read next line
        AK = fscanf(fid_AK,'%g',[C.retrieve_T(1)*C.retrieve_T(3)+C.retrieve_Q(1)*C.retrieve_Q(3)+C.retrieve_LWP,...
            C.retrieve_T(1)*C.retrieve_T(3)+C.retrieve_Q(1)*C.retrieve_Q(3)+C.retrieve_LWP]);
        nb_prof=nb_prof+1;
        
        R(nb_prof).nobs       = nb_prof;
        R(nb_prof).AK = AK;
    line = fgetl(fid_AK);
    if feof(fid_AK); break; end;
    line = fgetl(fid_AK);
    if feof(fid_AK); break; end;
    line = fgetl(fid_AK);
    if feof(fid_AK); break; end;         
    end
end

fclose(fid_AK);

end