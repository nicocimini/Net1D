function Qa = qspec_to_qabsolue(P,T,Qs)

% Function to convert specific humidity [kg/kg] into absolute humidity [kg/m3] 
% P pressure [Pa]
% T temperature in K
% Qs humidity in kg/kg

Ma=28.96;
Mv=18.015;
Ra=287.058;

Qa = (Qs .* P./T) ./ ( Ra * ( 1 - Qs*(1-Ma/Mv) ) );

return