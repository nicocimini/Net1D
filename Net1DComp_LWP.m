function [V,L] = Net1DComp_LWP(Tb23,Tb31)

% Pauline comment : If I understand well 
% outputs V = IWV in kg/m2
% outputs L= LWP in kg/m2
% inputs = brightness temperature measurements in K at 23 and 31 GHz

% DEFINE CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHAN #1
Tb01 = 8;  % offset
SV1 = 1.3; % V sensitivity
SL1 = 35;  % L sensitivity
% CHAN #2
Tb02 = 10; % offset
SV2 = 0.4; % V sensitivity
SL2 = 57;  % L sensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coefficients for Dual-channel retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
det = 1 / (SL2*SV1 - SL1*SV2);
c0 = (SL1*Tb02 - SL2*Tb01) * det;
c1 = SL2 * det;
c2 = -SL1 * det;
d0 = (SV2*Tb01 - SV1*Tb02) * det;
d1 = -SV2 * det;
d2 = SV1 * det;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computing Dual-channel retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = c0 + c1*Tb23 + c2*Tb31;
L = d0 + d1*Tb23 + d2*Tb31;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if L<0
    L=0;
end

return