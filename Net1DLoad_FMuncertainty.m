% Network 1DVAR+RTTOV retrieval: Load forward model uncertainty and sensitivity
%
% Net1DLoad_FMuncertainty load forward model uncertainty and sensitivity for 
% computing forward model error due to forward model uncertainty 

% Yet to be developed
% NB: need to prepare Kb and Sb (see ProfilingUncertainty_20to60GHz_V2_20170619.docx)
%
% Computing Model Error covariance matrix (YET TO BE DONE! Need to estimate Kb)
% Kb is the Jacobian of the measurement with respect to model parameter b) 
% i.e. Kb = dF(xa,b)/db(ba) - See GAIA_CLIM_GA2_WG2_T2.1.2_V2.pptx for reference 
% Sp = (Dy * Kb) * B * (Dy * Kb)t
% Sp = (Dy * Kb) * B * (Dy * Kb)';
    

function [Kb,Sb] = Net1DLoad_FMuncertainty(atm)

% qua dobbiamo semplicemente caricare dei file in cui abbiamo Kb (dipendente 
% dall'atmosfera) e Sb (fissi)
% I files vanno preparati in base alle stime di sensitivity to model
% parameters (vedi Sensitivity_20to60GHz_V9_20170614.docx)

% Kb deve avere dimensioni nchannel*nparameter e unit? K/dp dove dp indica la perturbazione ai vari parameter 
% Sb deve avere dimensioni nparameter*nparameter e unit? dp^2

return