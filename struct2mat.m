% This function tranforms an MxNx...-dimension field of an S-dimension 
% structure in a MxNx...xS matrix. 
%
% Es: 
%    MTX = struct2mat(S,field);
%
% Input: 
%       S: S-dimension structure
%   field: field of S, MxNx..-dimension matrix
%
% Output:
%     MTX: MxNx...xS matrix 
%
% Hystory:
%    2001/01: First created, Nico
%    2014/05: Completely rewritten, Nico
 
function [MTX] = struct2mat(S,field)

Ssize = length(S);
Fsize = size(getfield(S,{1},field));
Ndim = length(Fsize);

if max(Fsize)==1
   dimstr = '';
%   MTX = nan(1,Ssize);
else
   dimstr = repmat(':,',1,Ndim);
%   MTX = nan([Fsize Ssize]);
end

for i = 1:Ssize
   
   TMP = getfield(S,{i},field); 
   eval(['MTX(' dimstr num2str(i) ') = TMP;']);

end

MTX = squeeze(MTX);

return

% Older version
% function [MTX]=struct2mat(S,field)
% 
% for i=1:length(S)
%    
%    TMP=getfield(S,{i},field);   
%    dim=size(TMP);
%    if dim(2)>dim(1)
%       MTX(:,i)=TMP';      
%    else
%       MTX(:,i)=TMP;
%    end
%    
% end
% 
% return