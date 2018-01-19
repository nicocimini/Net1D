% Network 1DVAR+RTTOV retrieval: Compute IWV from retrieved Q profile
%
% Net1DComp_IWV computes IWV from the retrieved Q profile
% The Ret_H_kgm3 and IWV fields are added to R.
% The HTotErr field is added to A.
% The HObsErr, HSysUnc, IWV_rnd, and IWV_sys fields are added to E.
%
% Example:
%     [R,E] = Net1DComp_IWV(C,R,A,E);

function [R,A,E] = Net1DComp_IWV(C,X,R,A,E);

% Converting seconds from midnight to seconds from 1/1/1970
%time_bnds = C.sampling;
%time_offset = ( datenum(C.day_one(1),C.day_one(2),C.day_one(3)) - datenum(1970,1,1) ) * 24 * 3600;

nret = length(R);
%time = zeros(nret,1); 
for ir = 1:nret
    ibkg = X.ObsMatch(R(ir).nobs_prf);
%    R(ir).time = Otime(R(ir).nobs_tbs) + time_offset;
%    if R(ir).nite <= C.MaxIterations % convergence was reached
       % From specific (kg/kg) to absolute (kg/m3) humidity
       R(ir).Bkg_H_kgm3 = qspec_to_qabsolue(R(ir).Bkg_P_hPa*100,R(ir).Ret_T_K,R(ir).Bkg_Q_kgkg);
       R(ir).Ret_H_kgm3 = qspec_to_qabsolue(R(ir).Bkg_P_hPa*100,R(ir).Ret_T_K,R(ir).Ret_Q_kgkg);
       A(ir).HTotErr    = qspec_to_qabsolue(R(ir).Bkg_P_hPa*100,R(ir).Ret_T_K,A(ir).QTotErr);
       E(ir).HObsErr    = qspec_to_qabsolue(R(ir).Bkg_P_hPa*100,R(ir).Ret_T_K,E(ir).QObsErr); 
       E(ir).HSysUnc    = qspec_to_qabsolue(R(ir).Bkg_P_hPa*100,R(ir).Ret_T_K,E(ir).QSysUnc);
       % Computing IWV and associated uncertainty
       % Vertical integral of H[kg/m3] to compute [kg/m2]==[mm] (flip upside-down to have integral positive sign)
       z_m = flipud( X.Z(:,ibkg) ); 
       Bkg_H_kgm3 = flipud(R(ir).Bkg_H_kgm3);
       Ret_H_kgm3 = flipud(R(ir).Ret_H_kgm3);
       HObsErr = flipud(E(ir).HObsErr);
       HSysUnc = flipud(E(ir).HSysUnc);
       R(ir).Bkg_IWV_kgm2 = trapz(z_m, Bkg_H_kgm3 );     % IWV background value
       R(ir).Ret_IWV_kgm2 = trapz(z_m, Ret_H_kgm3 );     % IWV retrieved value
       E(ir).IWV_rnd = sqrt( trapz(z_m.^2,HObsErr.^2) ); % IWV rand uncertainty (based on ObsErr only)
       E(ir).IWV_sys = trapz(z_m,HSysUnc);               % IWV syst uncertainty (summed linearly)
%    end
end

return
