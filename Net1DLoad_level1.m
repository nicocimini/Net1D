% Network 1DVAR+RTTOV retrieval: Load level1 data
%
% Net1DLoad_level1 loads level 1 according to Config C structure, i.e. for
% different station, instrument, day,...

function [O,C] = Net1DLoad_level1(C);

% Initialization
O.y = []; % FixMe: null output
O.time = []; % FixMe: null output

% Searching for data files
lv1path = [C.datapath 'level1/' C.station_id '/'];
lv1head = ['*ps_' C.station(1:3) '_mwr00_l1_tb_v0*_'];
lv1date = [num2str(C.day_one(1)) twodigstr(C.day_one(2)) twodigstr(C.day_one(3))];
listofiles = dir([lv1path lv1head lv1date '*']);
lv1head_BL = ['*ps_' C.station(1:3) '_mwrBL00_l1_tb_v0*_'];
listofiles_BL = dir([lv1path lv1head_BL lv1date '*']);

% If file do not exist, return
if isempty(listofiles) & isempty(listofiles_BL);
    O.content = ['No lv1']; % Exit with null output
    return
end

% Loading is instrument independent (the output is not)
if ~isempty(listofiles_BL);
   % If BL file exist, take that one
   [O,C] = Net1DLoad_level1_BL(C,[lv1path listofiles_BL.name]);
   O.content = [C.instrument ' BL'];
   if C.std31fromZH % it ALSO load the Zenith obs and adopt the cloud flag (which should be more relialable)
      OZ = Net1DLoad_level1_ZH(C,[lv1path listofiles.name]);
      O.std31 = OZ.std31;
      O.cld31 = OZ.cld31;
      clear OZ;
   end
else
   % Otherwise zenith file
   [O,C] = Net1DLoad_level1_ZH(C,[lv1path listofiles.name]);
   O.content = [C.instrument ' ZH'];
   % NB: Here we need a check of C.elev_angles_degrees and C.elev_angles_how_many,
   % NB: as ZH files are not suited for retrievals with C.elev_angles_how_many > 1 or C.elev_angles_degrees ~= 90
   if C.elev_angles_how_many > 1 | C.elev_angles_degrees ~= 90
      fprintf(1,'Sorry, BL file is not available for this date.\n'); 
      fprintf(1,'Forcing C.elev_angles_how_many = 1 and C.elev_angles_degrees = 90 \n'); 
      fprintf(1,'Retrievals are computed with zenith only observations.\n');
      C.elev_angles_how_many = 1;
      C.elev_angles_degrees = 90;
   end
end        
    
end


