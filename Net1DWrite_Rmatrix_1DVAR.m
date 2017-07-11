function Net1DWrite_Rmatrix_1DVAR(C,Rmatrix)

%inputs for this function:
% C= config structure to run the 1DVAR
% ODVARpath = main path of the 1DVAR software
% R= the Rmatrix for the specified observation

if strcmp(C.instrument,'HATPRO')
    fid=fopen([C.ODVARpath,'HATPRO_COEFFS_DIR/','Rmatrix'],'w');
    fprintf(fid,'%s \n','HATPRO'); 
elseif strcmp(C.instrument,'MP3000A')
    fid=fopen([C.ODVARpath,'MP3000A_COEFFS_DIR/','Rmatrix'],'w');
    fprintf(fid,'%s \n','MP3000A'); 
end

if C.Rdiagonal
    % diagonal Rmatrix
    matrix_type=2;
else
    %Full Rmatrix
    matrix_type=1;
end

% total number of channels/TB observations that will be used in the 1DVAR
nb_chan=0;
for ang=1:C.elev_angles_how_many
    nb_chan=nb_chan+length(C.elev_angles_channum{ang});
end

if C.Rdiagonal
    % diagonal Rmatrix
    element_number=1; % only one band of diagonal
else
    %Full Rmatrix
    element_number=nb_chan;
end

R_inverse=0 ; % 0=R stored as inverse matrix, 1 = stored as direct matrix

fprintf(fid,'%d %d %d %d \n',matrix_type,nb_chan,element_number,R_inverse) ;

% channel number
for ang=1:C.elev_angles_how_many
     fprintf(fid,'%d ',C.elev_angles_channum{ang}(:)); 
end

fprintf(fid,'%s \n','  '); 

%computation of a full Rmatrix with all elevation angles
Rfull=zeros(nb_chan,nb_chan); % filled with zero by default
first_chan=1;
for ang=1:C.elev_angles_how_many
    last_chan=first_chan+length(C.elev_angles_channum{ang})-1;
    Rfull(first_chan:last_chan,first_chan:last_chan)=Rmatrix(C.elev_angles_channum{ang},C.elev_angles_channum{ang});
    first_chan=last_chan+1;
end

if C.Rdiagonal
    % diagonal R matrix only diagonal terms are extracted
    diag_R=diag(Rmatrix);
    for ang=1:C.elev_angles_how_many
        fprintf(fid,'%.2f ',diag_R(C.elev_angles_channum{ang}(:)));
    end
else
    % full R matrix we fill the file with all the values line by line 
    for line=1:nb_chan
        fprintf(fid,'%.2f ',Rfull(line,:)); 
    end
    
end

fclose(fid);
    
    

