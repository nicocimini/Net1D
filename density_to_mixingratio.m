% function to density of different species into
% mixing ratio in kg/kg
% inputs: rhox= density of species x in kg/m3
%         T= temperature in K
%         P= total pressure in Pa
%         rhov = water vapor density in kg/m3
% output= MR in kg/kg


function MR=density_to_mixingratio(rhox,T,P,rhov)

R=8.314510; %%% J/mol/K
rm=0.0289644 ; %%masse molaire de lair sec en kg/mol
Mv=18.015*1e-3;

% computation of the partial pressure of dry air by computing first
% partial pressure of water vapor

e=rhov*R.*T/Mv; % partial pressure of water vapor 
Pd=P-e; % partial pressure of dry air

MR=rhox*R.*T./(Pd*rm);
