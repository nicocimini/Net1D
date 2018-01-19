% Network 1DVAR+RTTOV retrieval: Compute retrieval error budget
%
% Net1DComp_errorbudget compute the error budget of the retrievals 
%
% NB: See GAIA_CLIM_GA2_WG2_T2.1.2_V2.pptx for reference 
% NB: I verified that A(ip).TTotErr = sqrt(diag(E(ip).SeT+E(ip).SaT) within 1e-4 K and 1e-6 kg/kg
% NB: The model error is still to be evaluated (need to evaluate Kb, i.e. Jacobian of the measurement with respect to model parameter b) 

function E = Net1DComp_errorbudget(C,R,B,AK,J,ybias)

np = length(AK);
nl = length(AK(1).AK);

if C.retrieve_LWP
    I = eye(nl-1);
else
    I = eye(nl);
end
Rinv = inv(R);
Binv = inv(B);

for ip = 1:np

    if C.retrieve_LWP
        K = J(ip).Jac(:,1:end-1); % nf * (nret*nl)
    else
        K = J(ip).Jac;
    end
    Kt = K';

    % Compute Dy (also called G) matrix
    % Dy = (B^-1 + Kt * R^-1 * K)^-1 * Kt * R^-1
    Dy = inv(Binv + Kt * Rinv * K) * Kt * Rinv;
    
    % Computing Observation Error covariance matrix
    % Se = Dy * R * Dyt
    Se = Dy * R * Dy';

    % Computing Smoothing Error covariance matrix
    % Ss = (AK - I) * B * (AK - I)t
    if ~C.retrieve_LWP
     Ss = (AK(ip).AK - I) * B * (AK(ip).AK - I)';
    else
     Ss = (AK(ip).AK(1:end-1,1:end-1) - I) * B * (AK(ip).AK(1:end-1,1:end-1) - I)';   
    end
    
    % Computing Model Error covariance matrix (YET TO BE DONE! Need to estimate Kb)
    % Kb is the Jacobian of the measurement with respect to model parameter b) 
    % i.e. Kb = dF(xa,b)/db(ba) - See GAIA_CLIM_GA2_WG2_T2.1.2_V2.pptx for reference 
    % [Kb,Sb] = Net1DLoad_FMuncertainty(atm);
    % P = (Dy * Kb) * Sb * (Dy * Kb)t
    % P = (Dy * Kb) * Sb * (Dy * Kb)';
    % NB: each column of (Dy * Kb) * Sb should contain the error contribution of each parameter (I'm guessing)

    % Store error covariance matrix and profiles
    E(ip).Sa = B; % a priori
    E(ip).Se = Se;
    E(ip).Ss = Ss;
    E(ip).BkgErr = sqrt( diag(B) );
    E(ip).ObsErr = sqrt( diag(Se) );
    E(ip).SmtErr = sqrt( diag(Ss) );
        
    % Systematic uncertainty
    E(ip).SysUnc = Dy * ybias;
    
    % Divide into T and Q
    iT = C.retrieve_T(3);
    iQ = C.retrieve_Q(3);
    iret = C.retrieve_T(1) + C.retrieve_Q(1);
    if iret == 2
       E(ip).SaT = B(1:iT,1:iT);
       E(ip).SeT = Se(1:iT,1:iT);
       E(ip).SsT = Ss(1:iT,1:iT);
       E(ip).TBkgErr = sqrt( diag(E(ip).SaT) );
       E(ip).TObsErr = sqrt( diag(E(ip).SeT) );
       E(ip).TSmtErr = sqrt( diag(E(ip).SsT) );
       E(ip).TSysUnc = E(ip).SysUnc(1:iT);       
       E(ip).SaQ = B(iT+1:iT+iQ,iT+1:iT+iQ);
       E(ip).SeQ = Se(iT+1:iT+iQ,iT+1:iT+iQ);
       E(ip).SsQ = Ss(iT+1:iT+iQ,iT+1:iT+iQ);
       E(ip).QBkgErr = sqrt( diag(E(ip).SaQ) );
       E(ip).QObsErr = sqrt( diag(E(ip).SeQ) );
       E(ip).QSmtErr = sqrt( diag(E(ip).SsQ) );       
       E(ip).QSysUnc = E(ip).SysUnc(iT+1:iT+iQ);       
    else
       if C.retrieve_T(1); 
           E(ip).SaT = B(1:iT,1:iT); 
           E(ip).SeT = Se(1:iT,1:iT); 
           E(ip).SsT = Ss(1:iT,1:iT); 
           E(ip).TBkgErr = sqrt( diag(E(ip).SaT) );
           E(ip).TObsErr = sqrt( diag(E(ip).SeT) );
           E(ip).TSmtErr = sqrt( diag(E(ip).SsT) );
           E(ip).TSysUnc = E(ip).SysUnc(1:iT);       
       end
       if C.retrieve_Q(1); 
           E(ip).SaQ = B(1:iQ,1:iQ); 
           E(ip).SeQ = Se(1:iQ,1:iQ); 
           E(ip).SsQ = Ss(1:iQ,1:iQ); 
           E(ip).QBckErr = sqrt( diag(E(ip).SaQ) );
           E(ip).QObsErr = sqrt( diag(E(ip).SeQ) );
           E(ip).QSmtErr = sqrt( diag(E(ip).SsQ) );
           E(ip).QSysUnc = E(ip).SysUnc(1:iQ);       
       end            
    end
    
end
    
return

