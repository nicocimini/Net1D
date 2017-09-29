% function to convert mixing ratio of different species into density in
% kg/m3
% inputs: MR= mixing ratio in kg/kg
%         T= temperature in K
%         P= total pressure in Pa
%         rhov = water vapor density in kg/m3 
% output: density in kg/m3

function density=mixingratio_to_density(MR,T,P,rhov)

R=8.314510; %%% J/mol/K
rm=0.0289644 ; %%masse molaire de lair sec en kg/mol
Mv=18.015*1e-3;

% computation of the partial pressure of dry air by computing first
% partial pressure of water vapor

e=rhov*R.*T/Mv; % partial pressure of water vapor 
Pd=P-e; % partial pressure of dry air

density=MR.*Pd*rm./(R*T) ;