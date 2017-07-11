function Net1DWrite_Channelchoice_1DVAR(C)

if strcmp(C.instrument,'HATPRO')
    fid=fopen([C.ODVARpath,'HATPRO_COEFFS_DIR/','ChannelChoice.dat'],'w');
elseif strcmp(C.instrument,'MP3000A')
    fid=fopen([C.ODVARpath,'MP3000A_COEFFS_DIR/','ChannelChoice.dat'],'w');
end

%total number of observations that will be processed
nb_chan=0;
for ang=1:C.elev_angles_how_many
    nb_chan=nb_chan+length(C.elev_angles_channum{ang});
end

fprintf(fid,'%d \n',nb_chan);

for obs=1:nb_chan
  fprintf(fid,'%d %d %d \n',obs,1023,3);  
end

fclose(fid);
