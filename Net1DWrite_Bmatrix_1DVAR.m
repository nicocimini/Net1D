function Net1DWrite_Bmatrix_1DVAR(X,C)

% function to write the background error covariance matrix to the file format asked by the 1DVAR
% inputs: the X structure for the background where we will extract the Bmatrix
% C= configuration file init strucure for the 1DVAR reprocessing

if strcmp(C.instrument,'HATPRO')
    fid=fopen([C.ODVARpath,'HATPRO_COEFFS_DIR/','Bmatrix'],'w');
elseif strcmp(C.instrument,'MP3000A')  
    fid=fopen([C.ODVARpath,'MP3000A_COEFFS_DIR/','Bmatrix'],'w');
end

fprintf(fid,'%s \n',' Background error covariance matrix');
fprintf(fid,'%s \n',' Background error covariance matrix');
fprintf(fid,'%d \n',length(X.B));

% full B matrix we fill the file with all the values line by line 
for line=1:length(X.B)
    fprintf(fid,'%e ',X.B(line,:)); 
    fprintf(fid,'%s \n','');
end

fprintf(fid,'%s \n',' Background error covariance matrix');
fprintf(fid,'%s \n',' Background error covariance matrix');
fprintf(fid,'%d \n',length(X.B));

% full B matrix we fill the file with all the values line by line 
for line=1:length(X.B)
    fprintf(fid,'%e ',X.B(line,:)); 
    fprintf(fid,'%s \n','');
end

fclose(fid);






