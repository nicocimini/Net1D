% Network 1DVAR+RTTOV retrieval: Load Jacobians
%
% Net1DLoad_1DVARout_Jacobians loads output 1DVAR data according to Config C 
% structure, i.e. for different station, day,...

function R = Net1DLoad_1DVARout_Jacobians(C,X)

% find output files
out_Jac = [C.ODVARpath_retrieval '/RetJacobian.out'];

% open output files
fid_Jac = fopen(out_Jac,'r');

nb_chan=0;
for ang=1:C.elev_angles_how_many
    nb_chan=nb_chan+length(C.elev_angles_channum{ang});
end

nb_prof=0;
while 1
    
    line = fgetl(fid_Jac);
  % disp(line)
    if strcmp(line(1:9),' ObNumber'); % here start a new retrieved profile
     %   nobs = str2num(line(22:27));
        line = fgetl(fid_Jac);
        %read next line
        Jac = fscanf(fid_Jac,'%g',[nb_chan,...
            C.retrieve_T(1)*C.retrieve_T(3)+C.retrieve_Q(1)*C.retrieve_Q(3)+C.retrieve_LWP]);
        nb_prof=nb_prof+1;
        
        R(nb_prof).nobs       = nb_prof;
        R(nb_prof).Jac = Jac;
    
        line = fgetl(fid_Jac);
        if feof(fid_Jac); break; end;
        line = fgetl(fid_Jac);
        if feof(fid_Jac); break; end;        
    end
end

fclose(fid_Jac);

end
