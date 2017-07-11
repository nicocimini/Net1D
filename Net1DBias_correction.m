% Network 1DVAR+RTTOV retrieval: Apply bias correction
%
% Net1DBias_correction apply bias correction according to Config C 
% structure, i.e. yes/no, different station,...

function O = Net1DBias_correction(C,O);

% Apply bias correction only if needed - otherwise exit with no action
if C.biascorrection

   % Load bias correction file
   load([C.biascorr_path C.biascorr_file]);

   % Apply bias correction
   bias_AROME_clearsky = bias_AROME_clearsky';
   O.y = O.y - bias_AROME_clearsky; % Tb_bc = Tb_o - bias = Tb_o - (Tb_o - Tb_b) = Tb_b

end

end
