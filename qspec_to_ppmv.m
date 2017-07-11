%%fonction pour convertir l'humidite specifique en ppmv

%%input: humidité spécifique en kg/kg
%%output: humidite en ppmv  

function qppmv=qspec_to_ppmv(q)

Ma=28.96;
Mv=18.015;

ratio=Mv/Ma;

v=q./(ratio*(1-q)+q);

%%% humid air molar mass

Mwet=(1-v)*Ma+v*Mv;

qppmv=q.*Mwet/Mv * 1.e6;