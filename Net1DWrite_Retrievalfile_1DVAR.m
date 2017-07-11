function Net1DWrite_Retrievalfile_1DVAR(C,nlev)

% function to write the retrieval file needed to configure the 1DVAR
% inputs: C structure from the configuration file
% nlev: number of vertical levels of the background profile

fid=fopen([C.ODVARpath,'Retrieval.NL'],'w');
fprintf(fid,'%s \n',''); 
fprintf(fid,'%s \n',''); 
fprintf(fid,'%s \n','! ******** THESE COMMENTS MAY NEED TO BE REMOVED IF COMPILING ***********'); 
fprintf(fid,'%s \n','! ******** WITH F90 (rather than F95) ***********************************'); 
fprintf(fid,'%s \n',''); 

%temperature retrieval
if C.retrieve_T(1)
    fprintf(fid,'%s %d %s %d %s %d \n','&Retrieval Temperature=', C.retrieve_T(2),',',C.retrieve_T(3),',',C.retrieve_T(2));
else
   fprintf(fid,'%s \n','&Retrieval Temperature= 0,0,0'); 
 
end

if C.retrieve_Q(1)
    fprintf(fid,'%s %d %s %d %s %d \n','Humidity=', C.retrieve_Q(2),',',C.retrieve_Q(3),',',C.retrieve_Q(2)+nlev);
else
   fprintf(fid,'%s \n','Humidity= 0,0,0'); 
end

   fprintf(fid,'%s \n','Surface_Temperature= 0,0 '); 
   fprintf(fid,'%s \n','Surface_Humidity= 0,0'); 
   fprintf(fid,'%s \n','Surface_Pressure= 0,0'); 
   fprintf(fid,'%s \n','Skin_Temperature= 0,0'); 
  
  if C.retrieve_LWP
     fprintf(fid,'%s \n',' Cloud_Liquid_Water= 1,0 ');  
  else
      fprintf(fid,'%s \n',' Cloud_Liquid_Water= 0,0 ');
  end
  
fprintf(fid,'%s \n','Integrated_Water_Vapor=0,0'); 
fprintf(fid,'%s \n','Cloud_Top_Pressure= 0,0'); 
fprintf(fid,'%s \n','Cloud_Fraction= 0,0'); 
fprintf(fid,'%s \n','Surface_UWind= 0,0 '); 
fprintf(fid,'%s \n','Surface_VWind= 0,0 /'); 

  
 

