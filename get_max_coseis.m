function nsta=get_max_coseis



%DMM 10/2011
%
%Simulate metrics to extract real-time coseismic offsets. This will fix
%offsets once the station witht he alrgest offset has finalized tot he
%offset.
%
%Requires you to save your MOVING AVERAGE filtered data to a structure
%called "coseis" witht he fields:
%
% coseis.N - North displacements
% coseis.E - East displacementsnt=get
% coseis.U - Up displacements


%Load data
%dir='/diego-local/Research/Data/Escenario2011/'
dir='/Users/dmelgarm/Research/Data/Tohoku/';
%dir='/diego-local/Research/Data/Toki/';
cd(dir)
%load d_MA50.mat
%load d_all_reduced_300.mat
%load toki300.mat
%load mx50.mat
%load tohoku_baseline_complete_waveforms.mat
load tohoku_gps_all.mat
%i=find(tohoku.complete==1);
% iname=tohoku.gname(i);
% for k=1:length(i)
%     ni(k)=find(coseis.gname==iname(k));
% end

% load tohokukal.mat
i=1:1:length(coseis.N)

%Run parameters
dt=0.01;   %sampling interval
Tavg=20;  %Average for final offset
thresh=0.0;   %minimum horizontal offset to be considered
p=0.1;     %Percent of peak variance for detection
Tvar=20/dt;   %Lagging variance window
Navg=floor(Tavg/dt);
%Limits for plotting later
xl=[0 500];
y1=[-5 5];
t1=50;
tmin=200

%coseis=tohoku;
% clear tohoku

% stdu=coseis.stdukal(i);
% stdn=coseis.stdnkal(i);
% stde=coseis.stdekal(i);
% lat=coseis.lat(i);
% lon=coseis.lon(i);
% n=coseis.nmakal(:,i)';
% e=coseis.emakal(:,i)';
% u=coseis.umakal(:,i)';
% T=(rnd(coseis.tk(:,i)*100)/100)';

% stdu=coseis.stdugps(n);
% stdn=coseis.stdngps(n);
% stde=coseis.stdegps(n);
% lat=coseis.lat(n);
% lon=coseis.lon(n);
% n=coseis.nmagps(:,n)';
% e=coseis.emagps(:,n)';
% u=coseis.umagps(:,n)';
% T=rnd(coseis.tgps(:,n))';

stdu=coseis.stdu(ni);
stdn=coseis.stdn(ni);
stde=coseis.stde(ni);
lat=coseis.lat(ni);
lon=coseis.lon(ni);
n=coseis.N(ni,:);
e=coseis.E(ni,:);
u=coseis.U(ni,:);
T=rnd(coseis.T(ni,:));

clear coseis
t=unique(T);
i=find(~isnan(t));
t=t(i)';
N=nan(size(T,1),size(T,2));
E=N;
U=N;
minT=min(T');
maxT=max(T');
nsta=size(E,1);
for k=1:nsta
   i1=find(t==minT(k));
   i2=find(t==maxT(k));
   j=find(~isnan(T(k,:)));
   j=j(end);
   N(k,i1:i2)=n(k,1:j);
   E(k,i1:i2)=e(k,1:j);
   U(k,i1:i2)=u(k,1:j);
end
T=repmat(t,nsta,1);
% N=coseis.N;
% E=coseis.E;
% U=coseis.U;
% T=coseis.T;
H=sqrt(E.^2+N.^2);
Hvar=movingSD(H,Tvar);   %Get trailing variance
mH=max(H');
i=find(isnan(H));
H(i)=0;
i=find(isnan(Hvar));
Hvar(i)=0;
iter=size(Hvar,2);
%get maximum var and horiz. displacement at each epoch
mH=zeros(size(Hvar));
mH(:,1)=H(:,1);
mHvar=zeros(size(Hvar));
mHvar(:,1)=Hvar(:,1);
for k=2:iter
    %Get maximum offsets
    mHvar_temp=Hvar(:,k);
    mH_temp=H(:,k);
    %Get maximum offsets
    gt=mH_temp>mH(:,k-1);  %If new epoch offset is alrger
    i=find(gt); %Which stations need updates
    mH(:,k)=mH(:,k-1); %these are the old maximum values
    mH(i,k)=mH_temp(i);  %And these are the one being updated
    %Get maximum variances
    gt=mHvar_temp>mHvar(:,k-1);  %If new epoch variance is alrger
    i=find(gt); %Which stations need updates
    mHvar(:,k)=mHvar(:,k-1); %these are the old maximum values
    mHvar(i,k)=mHvar_temp(i);  %And these are the one being updated
end
%Find when station with max offset drops below p% of Variance but
%respect horizontal offset threshold
detect=0;


[mH n_mH]=max(mH);
k=0;
while detect==0
    k=k+1;
    sta_mH=n_mH(k);  %Statio with maximum offset at that epoch
    if Hvar(Hvar(sta_mH,k)<p*mHvar(sta_mH,k))  & (t(k)>tmin)   %Final LARGEST offset reached
        %Compute offset over next Navg points
        Efinal(:,1)=mean(E(:,k:k+Navg),2);
        Nfinal(:,1)=mean(N(:,k:k+Navg),2);
        Ufinal(:,1)=mean(U(:,k:k+Navg),2);
        Tfinal(:,1)=T(:,k+Navg);
        detect=1;
        Hfinal=sqrt(Efinal.^2+Nfinal.^2);
        i=find(Hfinal<thresh); %Keep only stations voer threshold
        Efinal(i,:)=0;
        Nfinal(i,:)=0;
        Ufinal(i,:)=0;
        Hfinal(i,:)=0;
        %Get final variance function
        mHvar=Hvar(sta_mH,:);
        Tvar=T(sta_mH,:);
        finalvar=mHvar(k);
        finalTvar=Tfinal(1)-Tavg;
    end
    if k>iter
        detect=1;
        display('Reached end of run and could not compute final offsets!!!')
    end
end

 %Plot data
j=1:1:max(size(lat));
i=setxor(i,j); %Data over threshold

% figure
% subplot(2,1,1)
% Hmax=H(29,:);
% Tmax=T(29,:);
% Hmax=H(171,:);
% Tmax=T(171,:);
% Hmax=H(2,:);
% Tmax=T(2,:);
Hmax=H(1,:);
Tmax=T(1,:);
% 
% i=find(Tmax>=0 & Tmax<=275);
% Hmax=Hmax(i);
% Tmax=Tmax(i);
% i=find(Tvar>=0 & Tvar<=275);
% mHvar=mHvar(i);
% Tvar=Tvar(i);
ha = tight_subplot(4,1,0,0.1,0.1);

axes(ha(4));
[AX]=plotyy(Tmax,Hmax,Tvar,mHvar)
set(AX(1),'XTickLabel',[])
set(AX(2),'XTickLabel',[])
set(AX(1),'XLim',xl)
set(AX(2),'XLim',xl)
set(AX(1),'YLim',[0 6])
set(AX(2),'YLim',[0 0.2])
hold on
legend('Horizontal','Variance')
% stem(41,0.21,'--','LineWidth',2,'Color',[74/255 74/255 74/255])
% stem(49,0.21,'--','LineWidth',2,'Color',[74/255 74/255 74/255])
% stem(166,0.21,'--','LineWidth',2,'Color',[74/255 74/255 74/255])
% stem(186,0.21,'--','LineWidth',2,'Color',[74/255 74/255 74/255])
stem(t1,8,'--','LineWidth',2,'Color',[74/255 74/255 74/255])
stem(108,8,'--','LineWidth',2,'Color',[74/255 74/255 74/255])
stem(Tfinal(1)-Tavg,8,'--','LineWidth',2,'Color',[74/255 74/255 74/255])
stem(Tfinal(1),8,'--','LineWidth',2,'Color',[74/255 74/255 74/255])
grid on
xlabel('Seconds after origin time','FontSize',18)
xlim(xl)
set(get(AX(1),'Ylabel'),'String','Horizontal ([mm)','FontSize',18)
set(get(AX(2),'Ylabel'),'String','Variance (m^2)','FontSize',18)

axes(ha(1))
plot(T(i,:)',N(i,:)','LineWidth',2)
grid on
ylabel('North (m)')
hold on
plot(Tfinal',Nfinal','o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
plot(Tfinal'-Tavg,Nfinal','o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
plot([Tfinal (Tfinal-Tavg)]',[Nfinal Nfinal]','k','LineWidth',2)
plot([t1 t1],y1,'--k','LineWidth',2)
ylim([-3 0.5])
xlim([xl])
set(gca,'XTickLabel',[])
set(gca,'XTick',xl(1):20:xl(2))
axes(ha(2))
plot(T(i,:)',E(i,:)','LineWidth',2)
grid on
ylabel('East (m)')
hold on
plot(Tfinal',Efinal','o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
plot(Tfinal'-Tavg,Efinal','o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
plot([Tfinal (Tfinal-Tavg)]',[Efinal Efinal]','k','LineWidth',2)
plot([t1 t1],y1,'--k','LineWidth',2)
ylim([-0.5 6])
xlim([xl])
set(gca,'XTickLabel',[])
set(gca,'XTick',xl(1):20:xl(2))
axes(ha(3))
plot(T(i,:)',U(i,:)','LineWidth',2)
grid on
ylabel('Up (m)')
hold on
plot(Tfinal',Ufinal','o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
plot(Tfinal'-Tavg,Ufinal','o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
plot([Tfinal (Tfinal-Tavg)]',[Ufinal Ufinal]','k','LineWidth',2)
plot([t1 t1],y1,'--k','LineWidth',2)
ylim([-3 0.5])
xlim([xl])
set(gca,'XTickLabel',[])
set(gca,'XTick',xl(1):20:xl(2))


%Save
cd(dir)
coseis.N=Nfinal;
coseis.E=Efinal;
coseis.U=Ufinal;
coseis.T=Tfinal;
coseis.lat=lat;
coseis.lon=lon;
coseis.stdn=stdn;
coseis.stde=stde;
coseis.stdu=stdu;

jmax=find(Hfinal>thresh);
display([num2str(max(size(jmax))) ' stations over threshold.'])
pwd
save('tohoku_coseis.mat','coseis')
nsta=j;