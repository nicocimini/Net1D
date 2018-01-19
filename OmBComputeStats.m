% O minus B monitoring: compute O-B stats
%
% OmBComputeStats compute the O-B stats.
% It loads all the days from the selected period.
%
% Example:
%   OmBpath = '/Users/Nico/PROGETTI/TOPROF/SCIENCE/O-B/DATA/TOPROF/OB_bias/sir/';
%   datefrom = [2014 09 01];
%   dateto = [2014 09 09];
%   datebad = [2014 09 08; 2014 09 09];
%   [OmB,ALL] = OmBComputeStats(OmBpath,datefrom,dateto,datebad);
%   OmBPlotStats(OmB,OmBpath)

function [OmB,ALL] = OmBComputeStats(OmBpath,datefrom,dateto,datebad)

% If datebad is not provided, all are good
if nargin < 4
    datebad = [];
end

% Select days
Nstart = datenum(datefrom);
Nstop  = datenum(dateto);
Ndays = Nstart:Nstop;
if not(isempty(datebad))
   Nbad = datenum(datebad);
   Ndays(ismember(Ndays,Nbad)) = [];
end

% Load and concatenate
ALL = struct('time',[],'cld31',[],'std31',[],'TbObs',[],'TbBkg',[]);
for id = Ndays

    if exist([OmBpath 'OmB_' datestr(id,'yyyymmdd') '.mat'],'file');
       load([OmBpath 'OmB_' datestr(id,'yyyymmdd') '.mat']);
       
       %OmB.time(secs) -> computertime
       %thisdate = datevec(id);
       %julday   = date2jul(thisdate(1),thisdate(2),thisdate(3),0,0,OmB.time);
       time = ( OmB.time + (id - datenum(1970,1,1)) * 3600 * 24 ); % computertime
       ALL.time  = [ALL.time; time]; % Julian day
       ALL.cld31 = [ALL.cld31; OmB.cld31];
       ALL.std31 = [ALL.std31; OmB.std31];
       ALL.TbObs = cat(3,ALL.TbObs,OmB.TbObs);
       ALL.TbBkg = cat(3,ALL.TbBkg,OmB.TbBkg);
    
    end

end

ALL.station_id = C.station_id;
ALL.instrument = C.instrument;
ALL.channels   = OmB.channels;
ALL.angles_el  = OmB.angles_el;


% Compute Stats
clear OmB;
OmB.time = ALL.time;
diff = ALL.TbObs - ALL.TbBkg; 
OmB.diff = diff;
% select clear sky
cloudy = find(ALL.cld31);
diff(:,:,cloudy) = [];
OmB.mean = mean(diff,3); 
OmB.stdv = std(diff,0,3);
OmB.cld31= ALL.cld31;
OmB.std31= ALL.std31;
OmB.station_id = ALL.station_id;
OmB.instrument = ALL.instrument;
OmB.channels   = ALL.channels;
OmB.angles_el  = ALL.angles_el;


return