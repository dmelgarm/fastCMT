function prepare_elmayor(fsave,ma)

%First epoch is at 22:30:00 GPS time
%Origin time is 22:41:00 GPS time

gpsT=22*3600+30*60; %Time at start of record
originT=22*3600+41*60; %Origin time
originT=originT-gpsT; %relative origin time

cd('/diego-local/Research/Events/GPS')
load d.mat
%stafil=[13 27 33 46 60 61 62 63 64 65 66 67 69 70 71 72 77 81];
%stafil=1:1:104;
cd /diego-local/Research/Events/GPS
[lon lat a a a a]=textread('SEM_coseis_AZRYfix.txt','%f%f%f%f%f%s');
%Clean up
stafil=find(lon>-117.4 & lat<33.58);
lon=lon(stafil);
lat=lat(stafil);

E=E(:,stafil);
N=N(:,stafil);
U=U(:,stafil);
T=T(:,stafil);



%Clean up incomplete time series
[i c]=find(E==0);
c=unique(c);
c=setxor(c,1:1:max(size(stafil)));
E=E(:,c);
N=N(:,c);
U=U(:,c);
T=T(:,c);
lat=lat(c);
lon=lon(c);
%Get uncertainties
for k=1:min(size(T))
    stdn(k)=std(N(1:100,k)-mean(N(1:100,k)));
    stde(k)=std(E(1:100,k)-mean(E(1:100,k)));
    stdu(k)=std(U(1:100,k)-mean(U(1:100,k)));
end
%normstd=min([stdn stde stdu]);
normstd=1;
stdn=stdn/normstd;
stde=stde/normstd;
stdu=stdu/normstd;
%stdu=stdu*10000;
% stdn=ones(size(stdn));
% stde=ones(size(stdn));
% stdu=4*ones(size(stdn));
% stdn=stdn.^2;
% stde=stde.^2;
% stdu=stdu.^2;
% E=E';
% N=N';
% U=U';
% T=T';

iter=min(size(T));
T=T-originT;
for k=1:iter
    k
    i=find(T(:,k)>=0);
    nma=movavg(N(:,k),ma,ma);
    ema=movavg(E(:,k),ma,ma);
    uma=movavg(U(:,k),ma,ma);
    tma=T(:,k);
    Etemp(k,:)=ema(i);
    Ntemp(k,:)=nma(i);
    Utemp(k,:)=uma(i);
    Ttemp(k,:)=tma(i);
end
E=Etemp;
N=Ntemp;
U=Utemp;
T=Ttemp;

% mE=mean(E(:,1:500),2);
% mN=mean(N(:,1:500),2);
% mU=mean(U(:,1:500),2);
mE=E(:,1);
mN=N(:,1);
mU=U(:,1);
mE=repmat(mE,1,size(E,2));
mN=repmat(mN,1,size(N,2));
mU=repmat(mU,1,size(U,2));

E=E-mE;
N=N-mN;
U=U-mU;


%Filter by coseismic offset
H=sqrt(E.^2+N.^2);


%Keep only stuff after origin time
tfil=1:1:600;

coseis.E=E(:,tfil);
coseis.N=N(:,tfil);
coseis.U=U(:,tfil);
coseis.T=T(:,tfil);
coseis.lat=lat;
coseis.lon=lon;
coseis.stdn=stdn;
coseis.stde=stde;
coseis.stdu=stdu*1000;
cd('/diego-local/Research/Data/El Mayor')
save([fsave '.mat'],'coseis')


%plot'em
subplot(3,1,1)
plot(coseis.T',coseis.E','k','LineWidth',1.5)
grid on
ylabel('East (m)','FontSize',20)
subplot(3,1,2)
plot(coseis.T',coseis.N','k','LineWidth',1.5)
grid on
ylabel('North (m)','FontSize',20)
subplot(3,1,3)
plot(coseis.T',coseis.U','k','LineWidth',1.5)
grid on
ylabel('Up (m)','FontSize',20)
xlabel('Seconds after origin time','FontSize',20)