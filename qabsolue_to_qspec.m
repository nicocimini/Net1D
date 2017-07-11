function qspec=qabsolue_to_qspec(P,T,Q)

%%ùfonction pour convertir l'humidite absolue en kg/m3 en humidite specifique en kg/kg
%pression en Pa
%T en kg
% q en kg/m3

Ma=28.96;
Mv=18.015;
Ra=287.058;

C=Ra*Q.*T./P;

qspec=C./(1+C.*(1-Ma/Mv));