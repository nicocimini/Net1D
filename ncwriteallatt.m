% This function writes all the attributes of a given variable from a nc file into a new nc file.
%
% Usage:
% ncinf = ncinfo(ncinputfile);
% ncwriteallatt(ncoutputfile,ncinf.Variables(1));

function ncwriteallatt(ncoutputfile,variable)

natt = length(variable.Attributes);

for ia = 1:natt
 
    attname = variable.Attributes(ia).Name;
    if strcmp(attname,'_FillValue'); % for some reason the '_' at the beginning make Matlab to crash
        attname = attname(2:end);
    end
    ncwriteatt(ncoutputfile,variable.Name,attname,variable.Attributes(ia).Value);
    
end

return
