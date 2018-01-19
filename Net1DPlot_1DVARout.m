% Network 1DVAR+RTTOV retrieval: Plot 1DVAR output data
%
% Net1DLoad_plot1DVARout(C,X,R) plots output 1DVAR data according to Config C 
% structure, i.e. for different station, day,...

function Net1DPlot_1DVARout(C,O,X,R,E,A,AK,J)

% to be completed (and cleaned!)

% Figure output path
C.FIGSpath = [C.ODVARpath_output 'FIGS/' C.station_id '/'];

% this date
dateone = [num2str(C.day_one(1)) twodigstr(C.day_one(2)) twodigstr(C.day_one(3))]; 

% level 1 data
figure
plot(O.time/3600,squeeze(O.y(1,1,:))); grid on;
if C.biascorrection; head = 'Bias corrected '; else; head = ''; end;
xlabel('Time [h]'); ylabel([head 'Tb [K] @ ' num2str(O.channels(1)) ' GHz']); 
title([C.instrument ' (' C.station ') ' datestr(datenum(C.day_one))]);

% =======background temperature==========================================
figure;
set(gcf, 'PaperPositionMode', 'auto')
set(gcf, 'Renderer', 'zbuffer')
pcolor(X.time/3600,X.Z(:,1)/1e3,X.T); hc = colorbar;
set(get(hc,'Label'),'String','K');
set(gca,'FontSize',16)
caxis([270 280])
xlabel('Time [h]');
ylabel('Height [km asl]'); 
title(['AROME BCKG Ta (' C.station ') ' datestr(datenum(C.day_one))]);
%set(get(hc,'Label'),'String','K','FontWeight','bold');
shading flat; ylim([0 2]); xlim([0 24]); set(gca,'xtick',0:2:24);
if ~exist([C.FIGSpath 'T/BACKGROUND/'],'dir')
    mkdir([C.FIGSpath 'T/BACKGROUND/'])
end
saveas(gcf,[C.FIGSpath 'T/BACKGROUND/' C.station_id '_BACK_T_BL_' dateone],'png');

%==================background humidite=====================================
Ret_Q_kgkg = struct2mat(R,'Ret_Q_kgkg'); % read here to calculate the limits
figure 
set(gcf, 'Renderer', 'zbuffer')
pcolor(X.time/3600,X.Z(:,1)/1e3,X.Q); hc = colorbar;
set(get(hc,'Label'),'String','kg/kg');
xlabel('Time [h]'); ylabel('Height [km asl]'); 
set(gca,'FontSize',16)
caxis([min(min(min(X.Q)),min(min(Ret_Q_kgkg))) max(max(max(X.Q)),max(max(Ret_Q_kgkg)))])
title(['AROME BCKG Qa (' C.station ') ' datestr(datenum(C.day_one))]);
xlabel('Time [h]','FontSize',16); 
ylabel('Height [km asl]','FontSize',16); 
%set(get(hc,'Label'),'String','K','FontWeight','bold');
shading flat; ylim([0 10]); xlim([0 24]); set(gca,'xtick',0:2:24);
if ~exist([C.FIGSpath 'Q/BACKGROUND/'],'dir')
    mkdir([C.FIGSpath 'Q/BACKGROUND/'])
end
saveas(gcf,[C.FIGSpath 'Q/BACKGROUND/' C.station_id '_BACK_Q_' dateone],'png');


%=======Temperature retrieval=============================================
Ret_T_K = struct2mat(R,'Ret_T_K');
Back_T_K = struct2mat(R,'Bkg_T_K');
nobs = struct2mat(R,'nobs_prf');
nite = struct2mat(R,'nite');
% remove profiles with no convergence
noconv = find(nite>C.MaxIterations);
Ret_T_K(:,noconv) = [];
Back_T_K(:,noconv) = [];
nobs(noconv) = [];
% find gaps and put nans
gaps = find(diff(nobs)>1);
iadd = 0;
for ig = 1:length(gaps)
    ig1 = gaps(ig) + iadd;
%    nobs = [nobs(1:ig1) NaN nobs(ig1+1:end)];
    nobs = [nobs(1:ig1) nobs(ig1)+1 nobs(ig1+1:end)];
    Ret_T_K = [Ret_T_K(:,1:ig1) NaN(X.nlev,1) Ret_T_K(:,ig1+1:end)];
    Back_T_K = [Back_T_K(:,1:ig1) NaN(X.nlev,1) Back_T_K(:,ig1+1:end)];
    iadd = iadd + 1;
end

figure 
set(gcf, 'Renderer', 'zbuffer')
pcolor(O.time(nobs)/3600,X.Z(:,1)/1e3,Ret_T_K); hc = colorbar;
xlabel('Time [h]','FontSize',16); 
ylabel('Height [km asl]','FontSize',16); 
caxis([270 280])
set(get(hc,'Label'),'String','K','FontWeight','bold');
set(gca,'FontSize',16)
title(['1DVAR Ret Ta (' C.station ') ' datestr(datenum(C.day_one))]);
shading flat; ylim([0 2]); xlim([0 24]); set(gca,'xtick',0:2:24);
%format4paper(gcf);
if ~exist([C.FIGSpath 'T/1DVAR/'],'dir')
    mkdir([C.FIGSpath 'T/1DVAR/'])
end
saveas(gcf,[C.FIGSpath 'T/1DVAR/' C.station_id '_1DVAR_T_BL_' dateone],'png');

% Increment
figure 
set(gcf, 'Renderer', 'zbuffer')
pcolor(O.time(nobs)/3600,X.Z(:,1)/1e3,Ret_T_K-Back_T_K); hc = colorbar;
xlabel('Time [h]','FontSize',16); 
ylabel('Height [km asl]','FontSize',16); 
shading flat; ylim([0 2]);xlim([0 24]); set(gca,'xtick',0:2:24);
%caxis([270 280])
set(get(hc,'Label'),'String','K','FontWeight','bold');
set(gca,'FontSize',16)
title(['1DVAR Ret Ta - Back Ta (' C.station ') ' datestr(datenum(C.day_one))]);
if ~exist([C.FIGSpath 'T/1DVAR/'],'dir')
    mkdir([C.FIGSpath 'T/1DVAR/'])
end
saveas(gcf,[C.FIGSpath 'T/1DVAR/' C.station_id '_1DVAR_T_BL_Increment' dateone],'png');

%===humidity retrieval=====================================================
Back_Q_kgkg=struct2mat(R,'Bkg_Q_kgkg');
nobs = struct2mat(R,'nobs_prf');
nite = struct2mat(R,'nite');
% remove profiles with no convergence
noconv = find(nite>C.MaxIterations);
Ret_Q_kgkg(:,noconv) = [];
Back_Q_kgkg(:,noconv) = [];
nobs(noconv) = [];
% find gaps and put nans
gaps = find(diff(nobs)>1);
iadd = 0;
for ig = 1:length(gaps)
    ig1 = gaps(ig) + iadd;
%    nobs = [nobs(1:ig1) NaN nobs(ig1+1:end)];
    nobs = [nobs(1:ig1) nobs(ig1)+1 nobs(ig1+1:end)];
    Ret_Q_kgkg = [Ret_Q_kgkg(:,1:ig1) NaN(X.nlev,1) Ret_Q_kgkg(:,ig1+1:end)]; 
    Back_Q_kgkg = [Back_Q_kgkg(:,1:ig1) NaN(X.nlev,1) Back_Q_kgkg(:,ig1+1:end)];
    iadd = iadd + 1;
end

figure 
set(gcf, 'Renderer', 'zbuffer')
pcolor(O.time(nobs)/3600,X.Z(:,1)/1e3,Ret_Q_kgkg); hc = colorbar;
xlabel('Time [h]','FontSize',16); 
ylabel('Height [km asl]','FontSize',16); 
set(get(hc,'Label'),'String','kg/kg','FontWeight','bold');
set(gca,'FontSize',16)
title(['1DVAR Ret Qa (' C.station ') ' datestr(datenum(C.day_one))]);
shading flat;
ylim([0 10]); xlim([0 24]); set(gca,'xtick',0:2:24);
%caxis([0 8e-3])
caxis([min(min(min(X.Q)),min(min(Ret_Q_kgkg))) max(max(max(X.Q)),max(max(Ret_Q_kgkg)))])
%format4paper(gcf);
if ~exist([C.FIGSpath 'Q/1DVAR/'],'dir')
    mkdir([C.FIGSpath 'Q/1DVAR/'])
end
saveas(gcf,[C.FIGSpath 'Q/1DVAR/' C.station_id '_1DVAR_Q_' dateone],'png');

% increment
figure 
set(gcf, 'Renderer', 'zbuffer')
pcolor(O.time(nobs)/3600,X.Z(:,1)/1e3,Ret_Q_kgkg-Back_Q_kgkg); hc = colorbar;
xlabel('Time [h]','FontSize',16); 
ylabel('Height [km asl]','FontSize',16); 
shading flat; ylim([0 10]);xlim([0 24]); set(gca,'xtick',0:2:24);
%caxis([270 280])
set(get(hc,'Label'),'String','kg/kg','FontWeight','bold');
set(gca,'FontSize',16)
title(['1DVAR Ret Qa - Back Qa (' C.station ') ' datestr(datenum(C.day_one))]);
if ~exist([C.FIGSpath 'Q/1DVAR/'],'dir')
    mkdir([C.FIGSpath 'Q/1DVAR/'])
end
saveas(gcf,[C.FIGSpath 'Q/1DVAR/' C.station_id '_1DVAR_Q_Increment' dateone],'png');


%====LWP retrieval ========================================================
Ret_LWP=struct2mat(R,'Ret_LWP');
Back_LWP=struct2mat(R,'Bkg_LWP');
Ret_LWP(noconv) = [];
Back_LWP(noconv) = [];

iadd = 0;
for ig = 1:length(gaps)
    ig1 = gaps(ig) + iadd;
    Ret_LWP = [Ret_LWP(1:ig1) NaN Ret_LWP(ig1+1:end)];
    Back_LWP = [Back_LWP(1:ig1) NaN Back_LWP(ig1+1:end)];
    iadd = iadd + 1;
end

% Liquid water path
figure;
plot(O.time(nobs)/3600,Ret_LWP*1000,'x-k','Linewidth',2,'MarkerSize',10);hold on
plot(O.time(nobs)/3600,Back_LWP*1000,'x-r','Linewidth',2,'MarkerSize',10);hold on
%plot(X.time/3600,X.LWP*1000,'x-r','Linewidth',2,'MarkerSize',10);hold on
plot(O.time/3600,O.LWP*1000,'x-b','Linewidth',2,'MarkerSize',10);hold on
xlim([0 24]); set(gca,'xtick',0:2:24);
legend('1DVAR','Back.','MWR')
xlabel('Time [h]','FontSize',16); 
ylabel('Liquid water path [g/m2]','FontSize',16); 
title(['1DVAR LWP (' C.station ') ' datestr(datenum(C.day_one))],'FontSize',16);
set(gca,'FontSize',16)
if ~exist([C.FIGSpath 'LWP/'],'dir')
    mkdir([C.FIGSpath 'LWP/'])
end
saveas(gcf,[C.FIGSpath 'LWP/' C.station_id '_1DVAR_LWP_' dateone],'png');


%==== IWV retrieval =======================================================
Ret_IWV = struct2mat(R,'Ret_IWV_kgm2');
Bkq_IWV = struct2mat(R,'Bkg_IWV_kgm2');
Ret_IWV(noconv) = [];
Bkq_IWV(noconv) = [];

iadd = 0;
for ig = 1:length(gaps)
    ig1 = gaps(ig) + iadd;
    Ret_IWV = [Ret_IWV(1:ig1) NaN Ret_IWV(ig1+1:end)];
    Bkq_IWV = [Bkq_IWV(1:ig1) NaN Bkq_IWV(ig1+1:end)];
    iadd = iadd + 1;
end

% IWV
figure;
plot(O.time(nobs)/3600,Ret_IWV,'x-k','Linewidth',2,'MarkerSize',10);hold on
plot(O.time(nobs)/3600,Bkq_IWV,'x-r','Linewidth',2,'MarkerSize',10);hold on
xlim([0 24]); set(gca,'xtick',0:2:24);
legend('1DVAR','Bckgrd')
xlabel('Time [h]','FontSize',16); 
ylabel('IWV [kg/m^2]','FontSize',16); 
title(['1DVAR IWV (' C.station ') ' datestr(datenum(C.day_one))],'FontSize',16);
set(gca,'FontSize',16)
if ~exist([C.FIGSpath 'IWV/'],'dir')
    mkdir([C.FIGSpath 'IWV/'])
end
saveas(gcf,[C.FIGSpath 'IWV/' C.station_id '_1DVAR_IWV_' dateone],'png');


% ===estimated error (observation, smoothing, total, ...)=================
ip = 1;
figure
subplot(1,2,1)
plot(E(ip).TBkgErr,X.Z(:,ip)/1e3,'b','Linewidth',2); hold on
plot(E(ip).TObsErr,X.Z(:,ip)/1e3,'r','Linewidth',2); 
plot(E(ip).TSmtErr,X.Z(:,ip)/1e3,'m','Linewidth',2); 
%plot(sqrt(diag(X.B(1:60,1:60))),X.Z(:,ip)/1e3,'m','Linewidth',2); 
plot(A(ip).TTotErr,X.Z(:,ip)/1e3,'k','Linewidth',2); 
plot(E(ip).TSysUnc,X.Z(:,ip)/1e3,'k-.','Linewidth',2); 
% hold on; plot(sqrt(diag(E(ip).SeT+E(ip).SaT)),X.Z(:,ip)/1e3,'r--'); % to verify that A(ip).TTotErr = sqrt(diag(E(ip).SeT+E(ip).SaT)
xlabel('T [K]','FontSize',16); 
ylabel('Height [km asl]','FontSize',16); 
ylim([0 10]); grid on;
%set(gca,'xtick',[-0.2:0.2:1.0]);
set(gca,'FontSize',16)
legend('A priori','Obs','Smooth','Tot(rnd)','Tot(sys)');
subplot(1,2,2)
plot(E(ip).QBkgErr,X.Z(:,ip)/1e3,'b','Linewidth',2); hold on
plot(E(ip).QObsErr,X.Z(:,ip)/1e3,'r','Linewidth',2);
plot(E(ip).QSmtErr,X.Z(:,ip)/1e3,'m','Linewidth',2); 
plot(A(ip).QTotErr,X.Z(:,ip)/1e3,'k','Linewidth',2); 
plot(E(ip).QSysUnc,X.Z(:,ip)/1e3,'k-.','Linewidth',2);
% hold on; plot(sqrt(diag(E(ip).SeQ+E(ip).SaQ)),X.Z(:,ip)/1e3,'r--'); % to verify that A(ip).QTotErr = sqrt(diag(E(ip).SeQ+E(ip).SaQ)
xlabel('Q [kg/kg]','FontSize',16); ylim([0 10]); grid on;
%set(gca,'xtick',[0:2e-4:1e-3]);
set(gca,'FontSize',16)
legend('A priori','Obs','Smooth','Tot(rnd)','Tot(sys)');
if ~exist([C.FIGSpath 'ERROR/'],'dir')
    mkdir([C.FIGSpath 'ERROR/'])
end
saveas(gcf,[C.FIGSpath 'ERROR/1DVAR_Ret_err_prof_' dateone],'fig');
%format4paper(gcf);
saveas(gcf,[C.FIGSpath 'ERROR/1DVAR_Ret_err_prof_' dateone],'png');

% The same plot, but with absolute humidity (kg/m3)
% From specific (kg/kg) to absolute (kg/m3) humidity
Qa = qspec_to_qabsolue(R(ip).Bkg_P_hPa*100,R(ip).Ret_T_K,R(ip).Ret_Q_kgkg);
QBkgErr = qspec_to_qabsolue(R(ip).Bkg_P_hPa*100,R(ip).Ret_T_K,E(ip).QBkgErr);
QObsErr = qspec_to_qabsolue(R(ip).Bkg_P_hPa*100,R(ip).Ret_T_K,E(ip).QObsErr);
QSmtErr = qspec_to_qabsolue(R(ip).Bkg_P_hPa*100,R(ip).Ret_T_K,E(ip).QSmtErr);
QTotErr = qspec_to_qabsolue(R(ip).Bkg_P_hPa*100,R(ip).Ret_T_K,A(ip).QTotErr);
QSysUnc = qspec_to_qabsolue(R(ip).Bkg_P_hPa*100,R(ip).Ret_T_K,E(ip).QSysUnc);
figure
subplot(1,2,1)
plot(E(ip).TBkgErr,X.Z(:,ip)/1e3,'b','Linewidth',2); hold on
plot(E(ip).TObsErr,X.Z(:,ip)/1e3,'r','Linewidth',2); 
plot(E(ip).TSmtErr,X.Z(:,ip)/1e3,'m','Linewidth',2); 
plot(A(ip).TTotErr,X.Z(:,ip)/1e3,'k','Linewidth',2); 
plot(E(ip).TSysUnc,X.Z(:,ip)/1e3,'k-.','Linewidth',2); 
xlabel('T [K]','FontSize',16); 
ylabel('Height [km asl]','FontSize',16); 
ylim([0 10]); grid on;
%set(gca,'xtick',[0:0.2:1.0]);
set(gca,'FontSize',16)
legend('A priori','Obs','Smooth','Tot(rnd)','Tot(sys)');
subplot(1,2,2)
plot(QBkgErr,X.Z(:,ip)/1e3,'b','Linewidth',2); hold on
plot(QObsErr,X.Z(:,ip)/1e3,'r','Linewidth',2);
plot(QSmtErr,X.Z(:,ip)/1e3,'m','Linewidth',2); 
plot(QTotErr,X.Z(:,ip)/1e3,'k','Linewidth',2); 
plot(QSysUnc,X.Z(:,ip)/1e3,'k-.','Linewidth',2);
xlabel('Q [kg/m^3]','FontSize',16); ylim([0 10]); grid on;
%set(gca,'xtick',[0:2e-4:1e-3]);
set(gca,'FontSize',16)
legend('A priori','Obs','Smooth','Tot(rnd)','Tot(sys)');
if ~exist([C.FIGSpath 'ERROR/'],'dir')
    mkdir([C.FIGSpath 'ERROR/'])
end
saveas(gcf,[C.FIGSpath 'ERROR/1DVAR_Ret_err_prof_' dateone '_SpecHum'],'fig');
%format4paper(gcf);
saveas(gcf,[C.FIGSpath 'ERROR/1DVAR_Ret_err_prof_' dateone '_SpecHum'],'png');

% Absolute humidity only, but in (g/m3)
kgm3_to_gm3 = 1e3;
figure
plot(QBkgErr*kgm3_to_gm3,X.Z(:,ip)/1e3,'b','Linewidth',2); hold on
plot(QObsErr*kgm3_to_gm3,X.Z(:,ip)/1e3,'r','Linewidth',2);
plot(QSmtErr*kgm3_to_gm3,X.Z(:,ip)/1e3,'m','Linewidth',2); 
plot(QTotErr*kgm3_to_gm3,X.Z(:,ip)/1e3,'k','Linewidth',2); 
plot(QSysUnc*kgm3_to_gm3,X.Z(:,ip)/1e3,'k-.','Linewidth',2);
xlabel('Q [g/m^3]','FontSize',16); ylim([0 10]); grid on;
%set(gca,'xtick',[0:2e-4:1e-3]);
set(gca,'FontSize',16)
legend('A priori','Obs','Smooth','Tot(rnd)','Tot(sys)');


%==Profile with error bar================================================
ip = 1;
figure; % Temperature
%herrorbar(R(ip).Ret_T_K,X.Z(:,ip)/1e3,A(ip).TTotErr/2,'k'); hold on;
errorbar(R(ip).Ret_T_K,X.Z(:,ip)/1e3,A(ip).TTotErr/2,'k','horizontal'); hold on;
plot([R(ip).Ret_T_K-E(ip).TSysUnc R(ip).Ret_T_K+E(ip).TSysUnc],X.Z(:,ip)/1e3,'r--');
%plot(R(ip).Ret_T_K,X.Z(:,ip)/1e3,'k',[R(ip).Ret_T_K-A(ip).TTotErr/2 R(ip).Ret_T_K+A(ip).TTotErr/2],X.Z(:,ip)/1e3,'r:');
xlabel('T [K]','FontSize',16); 
ylabel('Height [km asl]','FontSize',16); 
ylim([0 10]); grid on;
set(gca,'FontSize',16);
legend('Rand','Syst');
title(['1DVAR Ret T (' C.station ') ' datestr( datenum([C.day_one])+O.time(R(1).nobs_tbs)/3600/24 )]);
saveas(gcf,[C.FIGSpath 'ERROR/1DVAR_Ret_err_prof_Prof_T' dateone],'png');

figure; % Specific humidity (kg/kg) 
%herrorbar(R(ip).Ret_Q_kgkg,X.Z(:,ip)/1e3,A(ip).QTotErr/2,'k'); hold on;
errorbar(R(ip).Ret_Q_kgkg,X.Z(:,ip)/1e3,A(ip).QTotErr/2,'k','horizontal'); hold on;
plot([R(ip).Ret_Q_kgkg-E(ip).QSysUnc R(ip).Ret_Q_kgkg+E(ip).QSysUnc],X.Z(:,ip)/1e3,'r--');
%plot(R(ip).Ret_Q_kgkg,X.Z(:,ip)/1e3,'k',[Ret_Q_kgkg(:,ip)-A(ip).QTotErr/2 R(ip).Ret_Q_kgkg+A(ip).QTotErr/2],X.Z(:,ip)/1e3,'r:');
xlabel('Q [kg/kg]','FontSize',16); ylim([0 10]); grid on;
ylabel('Height [km asl]','FontSize',16); 
ylim([0 10]); grid on;
set(gca,'FontSize',16)
legend('Rand','Syst');
title(['1DVAR Ret Q (' C.station ') ' datestr( datenum([C.day_one])+O.time(R(1).nobs_tbs)/3600/24 )]);
saveas(gcf,[C.FIGSpath 'ERROR/1DVAR_Ret_err_prof_Prof_Q' dateone],'png');

figure; % Absolute humidity (kg/m3) 
errorbar(Qa,X.Z(:,ip)/1e3,QTotErr,'k','horizontal'); hold on;
plot([Qa-QSysUnc Qa+QSysUnc],X.Z(:,ip)/1e3,'r--');
xlabel('Q [kg/m^3]','FontSize',16); ylim([0 10]); grid on;
ylabel('Height [km asl]','FontSize',16); 
ylim([0 10]); grid on;
set(gca,'FontSize',16)
legend('Rand','Syst');
title(['1DVAR Ret Abs Hum (' C.station ') ' datestr( datenum([C.day_one])+O.time(R(1).nobs_tbs)/3600/24 )]);
saveas(gcf,[C.FIGSpath 'ERROR/1DVAR_Ret_err_prof_Prof_AbsHum' dateone],'png');

       
%%%PLOT AVERAGING KERNEL
if C.retrieve_T
   figure;
   for lev = 1:C.retrieve_T(3)
       plot(AK(1).AK(lev,C.retrieve_T(2):C.retrieve_T(2)+C.retrieve_T(3)-1),X.Z(C.retrieve_T(2):C.retrieve_T(2)+C.retrieve_T(3)-1,1)/1e3,'k'); hold on
       % plot(AK(1).AK(C.retrieve_T(2):C.retrieve_T(2)+C.retrieve_T(3)-1,lev),X.Z(C.retrieve_T(2):C.retrieve_T(2)+C.retrieve_T(3)-1,1)/1e3,'k'); hold on
   end
       xlabel('T Averaging Kernel','FontSize',16);
    ylabel('Height [km asl]','FontSize',16);
    ylim([0 10])
    set(gca,'FontSize',16)
    title(['Averaging Kernel (' C.station ') ' datestr(datenum(C.day_one))]);
    if ~exist([C.FIGSpath 'T/AK/'],'dir')
        mkdir([C.FIGSpath 'T/AK/'])
    end
    saveas(gcf,[C.FIGSpath 'T/AK/' C.station_id '_AK_T_' dateone],'png');  
end

if C.retrieve_Q
   figure;
   for lev=1:C.retrieve_Q(3)
       plot(AK(1).AK(lev,C.retrieve_T(1)*C.retrieve_T(3)+1:C.retrieve_T(1)*C.retrieve_T(3)+C.retrieve_Q(3)),X.Z(C.retrieve_Q(2):C.retrieve_Q(2)+C.retrieve_Q(3)-1,1)/1e3,'k'); hold on
   end
     xlabel('Q Averaging Kernel','FontSize',16);
    ylabel('Height [km asl]','FontSize',16);
    ylim([0 10])
    set(gca,'FontSize',16)
    title(['Averaging Kernel (' C.station ') ' datestr(datenum(C.day_one))]);
    if ~exist([C.FIGSpath 'Q/AK/'],'dir')
        mkdir([C.FIGSpath 'Q/AK/'])
    end    
    saveas(gcf,[C.FIGSpath 'Q/AK/' C.station_id '_AK_Q_' dateone],'png'); 
end

    
%%%PLOT JACOBIANS
cmap=colormap(lines(C.channum));
if C.retrieve_T
   figure;
   for chan=1:C.channum
    plot(J(1).Jac(chan,1:C.retrieve_T(3)),X.Z(C.retrieve_T(2):C.retrieve_T(2)+C.retrieve_T(3)-1,1)/1e3,'Linewidth',2,'Color',cmap(chan,:)); hold on
   end
end    
    xlabel('T Jacobian','FontSize',16);
    ylabel('Height [km asl]','FontSize',16);
    h=legend(num2str(O.channels));
    legend boxoff
    set(h,'FontSize',12)
    ylim([0 10])
    set(gca,'FontSize',16)
    if ~exist([C.FIGSpath 'T/JAC/'],'dir')
        mkdir([C.FIGSpath 'T/JAC/'])
    end
    title(['Jacobian (' C.station ') ' datestr(datenum(C.day_one))]);
    saveas(gcf,[C.FIGSpath 'T/JAC/' C.station_id '_JAC_T_' dateone],'png');   

if C.retrieve_Q
   figure;
   for chan=1:C.channum
    %plot(AK(1).AK(lev,C.retrieve_T(1)*C.retrieve_T(3)+1:C.retrieve_T(1)*C.retrieve_T(3)+C.retrieve_Q(3)),X.Z(C.retrieve_Q(2):C.retrieve_Q(2)+C.retrieve_Q(3)-1,1)/1e3,'k','Linewidth',2); hold on
    plot(J(1).Jac(chan,C.retrieve_T(1)*C.retrieve_T(3)+1:C.retrieve_T(1)*C.retrieve_T(3)+C.retrieve_Q(3)),X.Z(C.retrieve_Q(2):C.retrieve_Q(2)+C.retrieve_Q(3)-1,1)/1e3,'k','Linewidth',2,'Color',cmap(chan,:)); hold on
   end
    xlabel('Q Jacobian','FontSize',16);
    ylabel('Height [km asl]','FontSize',16);
    h=legend(num2str(O.channels));
    legend boxoff
    set(h,'FontSize',12)    
    ylim([0 10])
    set(gca,'FontSize',16)
    title(['Jacobian (' C.station ') ' datestr(datenum(C.day_one))]);
    if ~exist([C.FIGSpath 'Q/JAC/'],'dir')
        mkdir([C.FIGSpath 'Q/JAC/'])
    end    
    saveas(gcf,[C.FIGSpath 'Q/JAC/' C.station_id '_JAC_Q_' dateone],'png');   
end

end




