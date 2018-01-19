% O minus B monitoring: Line up Obs and Bkg

function X = OmBLine_up_O_B(O,X)

% Settings
maxwndw = 900; % +/-15'

% Set dimensions
ntime = length(X.time);
nfrqs = length(O.channels);
nangs = length(O.angles_el);
X.channels = O.channels;
X.angles_el = O.angles_el;
X.TbObs = nan(nfrqs,nangs,ntime);
X.cld31 = nan(ntime,1);
X.std31 = nan(ntime,1);

for it = 1:ntime
    
    indx_O = find( abs(X.time(it) - O.time) < maxwndw ); 
    
    if isempty(indx_O)
        
        X.TbObs(:,:,it) = nan(nfrqs,nangs,1);
        X.cld31(it) = NaN;
        X.std31(it) = NaN;
        
    else
        %[it,indx_O]
        X.TbObs(:,:,it) = O.y(:,:,indx_O);
        X.cld31(it) = O.cld31(indx_O);
        X.std31(it) = O.std31(indx_O);
        
    end

end

% remove Nans (i.e. time slots where obs are not available)
indx = isnan(X.cld31);
X.time(indx) = [];
X.ps(indx) = [];
X.T(:,indx) = [];
X.P(:,indx) = [];
X.Z(:,indx) = [];
X.Q(:,indx) = [];
X.LWC(:,indx) = [];
X.LWP(indx) = [];
X.TbObs(:,:,indx) = [];
X.cld31(indx) = [];
X.std31(indx) = [];

return