%%% Mars 2016 P.Martinet : Creation
% This routine directly produces background file in the format expected by the NWPSAF 1D-Var
% input: P= Pressure profile in Pa
%        T= Temperature Profile in K
%        Q = Specific Humidity in kg/kg
%        Qsurf= Surface humidity in ppmv
%        Psurf = Surface pressur in Pa
%        CLW = cloud liquid water content in kg/kg
%        Fid = file identifier in which the data should be written
% number of the profile that is written among all the profiles that will be written in the file


function Net1DWrite_BACK_1DVAR(P,T,Q,CLW,Psurf,Qsurf,fid,nb_prof)

    fprintf(fid,'%s \n','x----------------------------------------------------------x ');
    fprintf(fid,'%s \n',sprintf('Profile #    %d Follows',nb_prof));
    fprintf(fid,'%s \n','x----------------------------------------------------------x ');
    for lev=1:length(P)
        fprintf(fid,'%f   %f   %e   %s   %e\n',P(lev)/100,T(lev),Q(lev),'0.00',CLW(lev));
    end
    fprintf(fid,'%s %f \n','Surface Temperature (K):        ',T(end));
    fprintf(fid,'%s %s \n','Surface Humidity (ppmv)  :   ', Qsurf);
    fprintf(fid,'%s %f \n','Skin Temperature (K):           ',T(end));
    fprintf(fid,'%s %f \n','Surface Pressure (hPa):         ',Psurf/100);
    fprintf(fid,'%s %s \n','10m U-Wind (m/s):                 ', '0.00');
    fprintf(fid,'%s %s \n','10m V-Wind (m/s):                 ', '0.00');
        
