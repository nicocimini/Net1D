% O minus B monitoring: Save Obs and Bkg output
%
% It stores the output OmB structure
%
% OmBSave_output(X,C)

function OmBSave_output(X,C)

out_path = [C.biascorr_path '/' C.station_id '/'];

fields = {'ps','T','P','Z','Q','LWC','LWP','B','nlev'};
OmB = rmfield(X,fields);

% naming output path and files
if ~exist(out_path,'dir'); mkdir(out_path); end;
out_file = ['OmB_' num2str(C.day_one(1)) twodigstr(C.day_one(2)) twodigstr(C.day_one(3))];
save([out_path out_file '.mat'],'C','OmB');


return