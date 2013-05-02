function [m lonmeca latmeca zmeca]=cmtgridmin(M,L,xyz)

% DMM 05/2010
%
% Find minimum moment tensor at each epoch make plots of its location
% (centroid) and write files for GMT movie script. Also make a misfit
% surface file for plotting in GMT and make displacements file.
%
% IN
% M - Momeont tensors
% L - Misfits
% xyz - Source point coordinates

% 
% OUT
% m - Moment tensor correpsonding to minimum misfit centroid
% lonmeca, latmeca, zmeca - Coordiantes of best centroid at each epoch

%Run parameters
dispflag=1;    %Make files for GMT plotting
misfitflag=0;  %Make files GMT plot of misfit surface
multiflag=1;   %Am I analyzing multiple point source data?


%dir='/diego-local/Research/Data/El Mayor/fastCMTelmay_L1R1';  %Where's the data from fastCMT?
dir='~/Research/Data/Tohoku/tohoku_RTO_kal';  %Where's the data from fastCMT?
%dir='/diego-local/Research/Data/Toki/fastCMTtoki_L1R2';  %Where's the data from fastCMT?
GMTpath='~/scripts/GMT/RTOkada/'  %Where is the GMT stuff?
t1=1;
%One inch = scale (m) on the plot
%scale=0.15; %Toki
%scale=0.04; %El Mayor  
scale=2; %Tohoku
%scale=3; %Mx
vscale=0.4;
%
%Filter by depth
% h1=0e3;
% h2=200e3;
% fildepth=find(xyz(:,3)>h1 & xyz(:,3)<h2);
% xyz=xyz(fildepth,:);
% M=M(:,:,:,fildepth);
% L=L(fildepth);



i=find(L==0);
L(i)=Inf;
[l n]=max(L);

% if size(L,1)>1
%     lon=xyz(n,1);
%     lat=xyz(n,2);
%     z=xyz(n,3)/1000;
% else
%     lon=xyz(1,1);
%     lat=xyz(1,2);
%     z=xyz(1,3)/1000;
%     n=ones(size(L));
% end
nt=size(L,2); %No of time slices
np=size(L,1); %No of source points% i=find(zgrid>20e3 & zgrid<60e3);
fildepth=1:1:np;

% zgrid=zgrid(i);
% latgrid=latgrid(i);
% longrid=longrid(i);
% i=find(latgrid>41 & latgrid<43);
% zgrid=zgrid(i);
% latgrid=latgrid(i);
% longrid=longrid(i);
% i=find(longrid>143 & longrid<145);
% zgrid=zgrid(i);
% latgrid=latgrid(i);
% longrid=longrid(i);


%Known solutions for plot

%Toki
% usgslat=linspace(42.171,42.171,nt);
% usgslon=linspace(143.719,143.719,nt);
% usgsz=linspace(13,13,nt);
% cmtlat=linspace(42.21,42.21,nt);
% cmtlon=linspace(143.84,143.84,nt);
% cmtz=linspace(28.2,28.2,nt);

%Tohoku
usgslat=linspace(37.92,37.92,nt);
usgslon=linspace(143.11,143.11,nt);
usgsz=linspace(19.5,19.5,nt);
cmtlat=linspace(37.52,37.52,nt);
cmtlon=linspace(143.05,143.05,nt);
cmtz=linspace(20.0,20.0,nt);

%El Mayor% i=find(zgrid>20e3 & zgrid<60e3);
% zgrid=zgrid(i);
% latgrid=latgrid(i);
% longrid=longrid(i);
% i=find(latgrid>41 & latgrid<43);
% zgrid=zgrid(i);
% latgrid=latgrid(i);
% longrid=longrid(i);
% i=find(longrid>143 & longrid<145);
% zgrid=zgrid(i);
% latgrid=latgrid(i);
% longrid=longrid(i);

% usgslat=linspace(32.192,32.192,nt);
% usgslon=linspace(-115.095,-115.095,nt);
% usgsz=linspace(15,15,nt);
% cmtlat=linspace(32.35,32.35,nt);
% cmtlon=linspace(-115.37,-115.37,nt);
% cmtz=linspace(12,12,nt);



% tplot=(1:1:nt);
% figure
% subplot(5,1,1)
% plot(tplot(t1:end),lon(t1:end),'LineWidth',2)
% grid on
% hold on
% plot(tplot(t1:end),usgslon(t1:end),'color',[160 32 40]/255,'LineWidth',2)
% plot(tplot(t1:end),cmtlon(t1:end),'color',[255 140 0]/255,'LineWidth',2)
% legend('GPS','USGS','Global CMT')
% ylabel('Lon (\circ)','FontSize',20)
% xlim([t1 300])
% subplot(5,1,2)
% plot(tplot(t1:end),lat(t1:end),'LineWidth',2)
% grid on
% hold on
% plot(tplot(t1:end),usgslat(t1:end),'color',[160 32 40]/255,'LineWidth',2)
% plot(tplot(t1:end),cmtlat(t1:end),'color',[255 140 0]/255,'LineWidth',2)
% ylabel('Lat (\circ)','FontSize',20)
% xlim([t1 300])
% subplot(5,1,3)
% plot(tplot(t1:end),z(t1:end),'LineWidth',2);
% grid on
% hold on
% plot(tplot(t1:end),usgsz(t1:end),'color',[160 32 40]/255,'LineWidth',2)
% plot(tplot(t1:end),cmtz(t1:end),'color',[255 140 0]/255,'LineWidth',2)
% ylabel('Depth (km)','FontSize',20)
% xlim([t1 300])
% 
% for k=1:nt
%     v(k)=L(n(k),k);
% end
% v(isinf(v))=0;

%Make misfit and moment tensor files (No displacements) for GMT movie
lon=xyz(:,1);
lat=xyz(:,2);
z=xyz(:,3);
if misfitflag==1
    lat2=linspace(min(lat),max(lat),200);
    lon2=linspace(min(lon),max(lon),200);
    L=L(1:end,:);
    [Lon Lat]=meshgrid(lon2,lat2);
    misfit=[];
    for k=1:nt
        k;
        v=L(:,k);
        V = griddata(lon,lat,v,Lon,Lat);
        [vmin row]=max(V);
        [vmin col]=max(vmin);
        row=row(col);
        glob_latmin(k)=Lat(row,col);
        glob_lonmin(k)=Lon(row,col);
        %Find closest one in real data
        [vmin i]=max(v);
        m(:,:,k)=M(:,:,k,i);
        depth(k)=z(i);
        misfit=[misfit lon lat v];
        
    end
    cd ~/scripts/GMT/fastCMT/Misfits
    save('misfit.xyz','misfit','-ascii')
end
%Make displacements fileif movieflaggeology quotes==1
%
%Data in variable d is ordered as follows
%1=station, 2=time, 3,4,5=observed data (Z,X,Y), 6,7,8=inverted
%data, 9,10, station coords
%x is north, y is east z is vertical with positive downwards
emisfit=[];
nmisfit=[];
umisfit=[];
umisfit_pos=[];
umisfit_neg=[];
vert=[];
%
if dispflag==1
    cd(dir);
    display('Creating flat files for movie...')
    k=1;
    for k=1:nt
        k
        v=L(:,k); %Load misfits for current time slice
        [vmin row]=max(v);   %Find position of min misfit
        load(['sdisp' num2str(fildepth(row)) '.mat'])
        load(['xyz' num2str(fildepth(row)) '.mat'])%Load displacements for min misfit position
        load(['moment' num2str(fildepth(row)) '.mat'])
        load(['mall' num2str(fildepth(row)) '.mat'])
        load(['source' num2str(fildepth(row)) '.mat'])
        %moment=moment(:,:,row);
        %Correct timesedit batch
        d(:,2)=d(:,2);
        lone=source(1);   %Source coords
        late=source(2);% stdu=coseis.stdukal(i);
% stdn=coseis.stdnkal(i);
% stde=coseis.stdekal(i);
% lat=coseis.lat(i);
% lon=coseis.lon(i);
% n=coseis.nmakal(:,i)';
% e=coseis.emakal(:,i)';
% u=coseis.umakal(:,i)';

        ze=source(3);
        lonmeca(k)=lone;
        latmeca(k)=late;
        zmeca(k)=ze;
        ts=unique(d(:,2));
        %Count max no. of stations
        nsta=size(unique(d(:,1)),1)-1;
        is=find(d(:,2)==k);
        ista=d(is,1);
        nsta=size(ista,1);
        if nsta>0 %if there is data to be saved
            lonsta=d(is,9);
            latsta=d(is,10);
            S(1:nsta,1+4*(k-1))=lonsta; %Station coords.
            S(1:nsta,2+4*(k-1))=latsta;
            vert(1:nsta,1+4*(k-1))=lonsta; %Station coords verticals.
            vert(1:nsta,2+4*(k-1))=latsta;
            vertS(1:nsta,1+4*(k-1))=lonsta+0.1; %Station coords verticals.
            vertS(1:nsta,2+4*(k-1))=latsta;
            Ss(1:nsta,1+4*(k-1))=S(1:nsta,1+4*(k-1));  %ditto...
            Ss(1:nsta,2+4*(k-1))=S(1:nsta,2+4*(k-1));
            S(1:nsta,3+4*(k-1))=rad2deg(atan2(d(is,4),d(is,5)));  %Coseismic azimuth (from East!!!)
            S(1:nsta,4+4*(k-1))=sqrt(d(is,5).^2+d(is,4).^2)/scale;  %Coseismic amplitude
            vert(1:nsta,3+4*(k-1))=90;
            vertS(1:nsta,3+4*(k-1))=90;
            Ss(1:nsta,3+4*(k-1))=rad2deg(atan2(d(is,7),d(is,8)));  %Coseismic azimuth, synthetic
            Ss(1:nsta,4+4*(k-1))=sqrt(d(is,8).^2+d(is,4).^2)/scale; %Coseismic amplitude, synthetic
            vert(1:nsta,4+4*(k-1))=-d(is,3)/vscale;  %Observed vertical
            vertS(1:nsta,4+4*(k-1))=-d(is,6)/vscale;  %Synthetic vertical
            %Fix verticals
            ineg=(vert(:,4)<0);
            vert(ineg,3+4*(k-1))=270;
            vert(ineg,4+4*(k-1))=-vert(ineg,4+4*(k-1));
            ineg=(vertS(:,4)<0);
            vertS(ineg,3+4*(k-1))=270;
            vertS(ineg,4+4*(k-1))=-vertS(ineg,4+4*(k-1));
            %Deltas are observed - synthetic
            deltau=-d(is,3)+d(is,6);
            deltan=d(is,4)-d(is,7);% stdu=coseis.stdukal(i);
% stdn=coseis.stdnkal(i);
% stde=coseis.stdekal(i);
% lat=coseis.lat(i);
% lon=coseis.lon(i);
% n=coseis.nmakal(:,i)';
% e=coseis.emakal(:,i)';
% u=coseis.umakal(:,i)';

            deltae=d(is,5)-d(is,8);
            %Recompute variances
            dobs=[d(is,3)' d(is,4)' d(is,5)'];
            dsyn=[d(is,6)' d(is,7)' d(is,8)'];
            it=find(dobs~=0);
            v2(k)=vmin;
%             if ~isempty(it)
%                 v2(k)=(1-((sum((dobs(it)-dsyn(it)).^2)/(sum(dobs(it).^2)))))*100;
%             else
%                 v2(k)=NaN;
%             end
            %Sort out positive and negative misfits
            mscale=0.05/0.5;  %0.05cm misfit is 0.5in in plot
            ip=find(deltau>=0);
            in=find(deltau<0);
            deltau_pos=deltau/mscale;
            deltau_neg=abs(deltau)/mscale;
            deltau_pos(in)=0;  %Set to zero to distinguish between symbol types later on
            deltau_neg(ip)=0;
        end    %Otherwise keep zeros
        %time vector
        ts(k)=k;
        if multiflag==1
            Mo=0;
            for km=1:size(M,3)
               Mo(km)=norm(M(:,:,km,row),'fro')/sqrt(2);
               Mcurrent=M(:,:,km,row);
            end
            Mw(k)=0.67*(log10(sum(Mo))-9.1);
        else    
            Mcurrent=M(:,:,k,row);  %Load current moment tensor
            m(:,:,k)=Mcurrent;
            Mo=norm(Mcurrent,'fro')/sqrt(2);  %Scalar momoent
            Mw(k)=0.67*(log10(Mo)-9.1);
        end
    end
    %Save displacements
    cd(GMTpath);
    save('vert.obs','vert','-ascii')
    save('vert.syn','vertS','-ascii')
    save('coseis.obs','S','-ascii');
    save('coseis.syn','Ss','-ascii');
    save('moment.xyz','moment','-ascii');
    save('grdsearch.xyz','xyz','-ascii')
    fid = fopen('coseis.time', 'w');
    fprintf(fid, '%3i\n', ts);
    fclose(fid);
    fid = fopen('Mw.dat', 'w');
    fprintf(fid, '%1.2f\n', Mw);
    fclose(fid);
    fid = fopen('VR.dat', 'w');
    fprintf(fid, '%4.2f\n', L);
    fclose(fid);
    if multiflag==1
        m=Mcurrent;
        nepi=max(size(xyz));
        nepi=floor(nepi/2);
        lonmeca=xyz(nepi,1);
        latmeca=xyz(nepi,2);
        zmeca=xyz(nepi,3);
        [lonmeca latmeca]=moment_centroid(moment);
        make_psmeca(m*1e7,lonmeca,latmeca,zmeca);
    else 
        make_psmeca(m*1e7,lonmeca,latmeca,zmeca);
    end
end
% subplot(5,1,5)
% plot(tplot,v2,'LineWidth',2)
% %ylim([0 max(v)*1.1])
% grid on
% ylabel('Variance Reduction (%)','FontSize',20)
% xlabel('Seconds after origin time','FontSize',20)
% xlim([t1 300])
% Mw(isinf(Mw))=0;
% subplot(5,1,4)
% plot(tplot,Mw,'LineWidth',2)
% ylim([5 10])
% xlim([t1 300])
% grid on
% ylabel('M_w','FontSize',20)

row
max(v)
display('Making report...')

  
fastCMTreport(m,lonmeca,latmeca,zmeca,L,139)
