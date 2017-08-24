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

% background temperature
figure;
set(gcf, 'PaperPositionMode', 'auto')
set(gcf, 'Renderer', 'zbuffer')
pcolor(X.time/3600,X.Z(:,1)/1e3,X.T); colorbar;
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

% background humidite
Ret_Q_kgkg = struct2mat(R,'Ret_Q_kgkg');
figure 
set(gcf, 'Renderer', 'zbuffer')
pcolor(X.time/3600,X.Z(:,1)/1e3,X.Q); colorbar;
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


% retrieval
Ret_T_K = struct2mat(R,'Ret_T_K');
nobs = struct2mat(R,'nobs_prf');
nite = struct2mat(R,'nite');
% remove profiles with no convergence
noconv = find(nite>C.MaxIterations);
Ret_T_K(:,noconv) = [];
nobs(noconv) = [];
% find gaps and put nans
gaps = find(diff(nobs)>1);
iadd = 0;
for ig = 1:length(gaps)
    ig1 = gaps(ig) + iadd;
%    nobs = [nobs(1:ig1) NaN nobs(ig1+1:end)];
    nobs = [nobs(1:ig1) nobs(ig1)+1 nobs(ig1+1:end)];
    Ret_T_K = [Ret_T_K(:,1:ig1) NaN(X.nlev,1) Ret_T_K(:,ig1+1:end)]; 
    iadd = iadd + 1;
end

figure 
set(gcf, 'Renderer', 'zbuffer')
pcolor(O.time(nobs)/3600,X.Z(:,1)/1e3,Ret_T_K); hc = colorbar;
xlabel('Time [h]','FontSize',16); 
ylabel('Height [km asl]','FontSize',16); 
caxis([270 280])
%set(get(hc,'Label'),'String','K','FontWeight','bold');
set(gca,'FontSize',16)
title(['1DVAR Ret Ta (' C.station ') ' datestr(datenum(C.day_one))]);
shading flat; ylim([0 2]); xlim([0 24]); set(gca,'xtick',0:2:24);
%format4paper(gcf);
if ~exist([C.FIGSpath 'T/1DVAR/'],'dir')
    mkdir([C.FIGSpath 'T/1DVAR/'])
end
saveas(gcf,[C.FIGSpath 'T/1DVAR/' C.station_id '_1DVAR_T_BL_' dateone],'png');

%humidite retrieval
nobs = struct2mat(R,'nobs_prf');
nite = struct2mat(R,'nite');
% remove profiles with no convergence
noconv = find(nite>C.MaxIterations);
Ret_Q_kgkg(:,noconv) = [];
nobs(noconv) = [];
% find gaps and put nans
gaps = find(diff(nobs)>1);
iadd = 0;
for ig = 1:length(gaps)
    ig1 = gaps(ig) + iadd;
%    nobs = [nobs(1:ig1) NaN nobs(ig1+1:end)];
    nobs = [nobs(1:ig1) nobs(ig1)+1 nobs(ig1+1:end)];
    Ret_Q_kgkg = [Ret_Q_kgkg(:,1:ig1) NaN(X.nlev,1) Ret_Q_kgkg(:,ig1+1:end)]; 
    iadd = iadd + 1;
end

figure 
set(gcf, 'Renderer', 'zbuffer')
pcolor(O.time(nobs)/3600,X.Z(:,1)/1e3,Ret_Q_kgkg); hc = colorbar;
xlabel('Time [h]','FontSize',16); 
ylabel('Height [km asl]','FontSize',16); 
%set(get(hc,'Label'),'String','K','FontWeight','bold');
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

Ret_LWP=struct2mat(R,'Ret_LWP');
Back_LWP=struct2mat(R,'Bkg_LWP');
Ret_LWP(noconv) = [];
Back_LWP(noconv) = [];

iadd = 0;
for ig = 1:length(gaps)
    ig1 = gaps(ig) + iadd;
%    nobs = [nobs(1:ig1) NaN nobs(ig1+1:end)];
 %   nobs = [nobs(1:ig1) nobs(ig1)+1 nobs(ig1+1:end)];
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

% estimated error (observation, smoothing, total, ...)
ip = 1;
figure
subplot(1,2,1)
plot(E(ip).TObsErr,X.Z(:,ip)/1e3,'r',E(ip).TSmtErr,X.Z(:,ip)/1e3,'b',A(ip).TTotErr,X.Z(:,ip)/1e3,'k',sqrt(diag(X.B(1:60,1:60))),X.Z(:,ip)/1e3,'m'); 
% hold on; plot(sqrt(diag(E(ip).SeT+E(ip).SaT)),X.Z(:,ip)/1e3,'r--'); % to verify that A(ip).TTotErr = sqrt(diag(E(ip).SeT+E(ip).SaT)
xlabel('T [K]','FontSize',16); 
ylabel('Height [km asl]','FontSize',16); 
ylim([0 10]); grid on;
set(gca,'xtick',[0:0.2:1.0]);
legend('Obs','Smooth','Tot','A priori');
subplot(1,2,2)
plot(E(ip).QObsErr,X.Z(:,ip)/1e3,'r',E(ip).QSmtErr,X.Z(:,ip)/1e3,'b',A(ip).QTotErr,X.Z(:,ip)/1e3,'k',sqrt(diag(X.B(61:120,61:120))),X.Z(:,ip)/1e3,'m'); 
% hold on; plot(sqrt(diag(E(ip).SeQ+E(ip).SaQ)),X.Z(:,ip)/1e3,'r--'); % to verify that A(ip).QTotErr = sqrt(diag(E(ip).SeQ+E(ip).SaQ)
xlabel('Q [kg/kg]'); ylim([0 10]); grid on;
set(gca,'xtick',[0:2e-4:1e-3]);
legend('Obs','Smooth','Tot','A priori');
if ~exist([C.FIGSpath 'ERROR/'],'dir')
    mkdir([C.FIGSpath 'ERROR/'])
end
saveas(gcf,[C.FIGSpath 'ERROR/1DVAR_Ret_err_prof_' dateone],'fig');
%format4paper(gcf);
saveas(gcf,[C.FIGSpath 'ERROR/1DVAR_Ret_err_prof_' dateone],'png');


% From Pauline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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




