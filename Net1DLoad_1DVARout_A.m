% Network 1DVAR+RTTOV retrieval: Load A matrix
%
% Net1DLoad_1DVARout_A loads output 1DVAR A matrix according to Config C 
% structure, i.e. for different station, day,...

function [A,Abkg] = Net1DLoad_1DVARout_A(C,X)

% find output files
out_A = [C.ODVARpath '/A-Matrix.out'];

% open output files
fid_A = fopen(out_A,'r');
  
npro = 0;
while 1
    
    line = fgetl(fid_A);
    if length(line) > 1

    if strcmp(line(1:12),' Calculated '); % the very first one is for background of first observation (just for example?)
        
        line = fgetl(fid_A);
        Abkg = fscanf(fid_A,'%g',[X.nlev,X.nlev]); % MTX = MTX'; % Need debuging? See the modification below by Pauline
        
    elseif strcmp(line(1:12),' Observation'); % here start a new retrieved profile
        
        nobs = str2num(line(30:34));
        line = fgetl(fid_A);
        %A1 = fscanf(fid_A,'%g',[X.nlev,X.nlev]); % MTX = MTX';
        A1 = fscanf(fid_A,'%g',[C.retrieve_T(1)*C.retrieve_T(3)+C.retrieve_Q(1)*C.retrieve_Q(3)+C.retrieve_LWP,...
            C.retrieve_T(1)*C.retrieve_T(3)+C.retrieve_Q(1)*C.retrieve_Q(3)+C.retrieve_LWP]);          
        npro = npro + 1;        
        A(npro).nobs = nobs;
        A(npro).A = A1;
        % Nico
        iT = C.retrieve_T(3);
        iQ = C.retrieve_Q(3);
        iret = C.retrieve_T(1) + C.retrieve_Q(1);
        if iret == 2
           A(npro).AT = A1(1:iT,1:iT);
           A(npro).AQ = A1(iT+1:iT+iQ,iT+1:iT+iQ);
           A(npro).TTotErr = sqrt( diag( A1(1:iT,1:iT) ) );
           A(npro).QTotErr = sqrt( diag( A1(iT+1:iT+iQ,iT+1:iT+iQ) ) );           
        else
            if C.retrieve_T(1); A(npro).AT = A1(1:iT,1:iT); A(npro).TTotErr = sqrt( diag( A1(1:iT,1:iT) ) ); end
            if C.retrieve_Q(1); A(npro).AQ = A1(1:iQ,1:iQ); A(npro).QTotErr = sqrt( diag( A1(1:iQ,1:iQ) ) ); end            
        end
        %
    end

    end
    
    if feof(fid_A); break; end;
    
end

fclose(fid_A);


return