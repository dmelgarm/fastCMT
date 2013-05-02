function batch_green(fgreen,data,fsave)

% DMM 08/2011
%
% Make Green functions for multiple inversion nodes
%
% IN:
% 
% fgreen ~ Name of Green function file computed in EDGRN
% data ~ Structure containing the data and metadata
% fsave ~ Filename for Green functions matrix (absolute path)

%Paths where EDGRN GFs are stored
pathgf='/diego-local/scripts/Fortran/ED/grnfcts'
%Path where data is stored
pathdata='/diego-local/Research/Data/Tohoku'
cd(pathdata)
%Compute GF derivatives?
derivflag=0;
dGF=2000;  %GF Spacing

% % Make inversion nodes
% %[longrid latgrid zgrid]=textread('fault.xyz','%f%f%f');
% %[longrid latgrid zgrid]=textread('slab1.0.xyz','%f%f%f');
% [longrid latgrid zgrid]=textread('slab1.0small.xyz','%f%f%f');
% longrid=longrid';
% latgrid=latgrid';
% zgrid=zgrid;
% zgrid=-zgrid';
% % i=find(longrid>134 & longrid<180);
% % zgrid=zgrid(i);
% % latgrid=latgrid(i);
% % longrid=longrid(i);
% % i=find(latgrid>40 & latgrid<44);
% % i=find(zgrid>150);
% % zgrid=zgrid(i);
% % latgrid=latgrid(i);
% % longrid=longrid(i);
% display('Creating grid on plate interface...')
% ________________
% 
% % % % Prism around epicentre
% % [longrid latgrid zgrid]=makegrid(-115.28,32.192,10e3,0.1,0.1,2e3,14,14,10);
% % %[longrid latgrid zgrid]=makegrid(144.00,42.2,35e3,0.2,0.2,3e3,14,14,10);
% % % [longrid latgrid zgrid]=makegrid(143.84,42.21,27e3,0.3,0.3,3e3,0,0,0)
% % [longrid latgrid zgrid]=makegrid(142.3,38.3,20e3,0.1,0.1,2e3,10,10,19);
% [longrid latgrid zgrid]=makegrid(142.4,38.3,30e3,0.1,0.1,2e3,20,20,30);
% % %[longrid latgrid zgrid]=makegrid(-100,16.5,20e3,0.2,0.2,2e3,16,16,16ongrid latgrid zgrid]=makegrid(37,138,50e3,0.01,0.01,1e3,1500,1500,100);
% i=find(zgrid>1e3);
% longrid=longrid(i);
% latgrid=latgrid(i);
% zgrid=zgrid(i);
% display('Creating prism grid around origin...')
% ____________

% Single sources
% % Toki
% longrid=[143.84]; 
% latgrid=[42.21];
% zgrid=[28e3];
% %Escenario2011
% latgrid= 16.5900;
% longrid= -99.0400;
% zgrid= 17.8200e3;
% % % El May
% % % longrid=[-115.37];
% % % latgrid=[32.35];
% % % zgrid=[8e3];
% Tohoku
longrid=142.369;
latgrid=38.32;
zgrid=20e3;
% display('No grid created, single point inversion...')
% _____________

% Multiple point sources
% longrid=[141.5 141.75 142.00 142.25 142.50 142.75 143.00 143.25 143.50 143.75 144.0]; 
% latgrid=[34.75 35.5 36.00 36.75 37.25 38.00 38.75 39.50 40.00 40.50 41.25];
% % zgrid=[20e3 20e3 20e3 20e3 20e3 20e3 20e3 20e3 20e3 20e3 20e3];
% N=50;
% longrid=linspace(141.5,144,N);
% latgrid=linspace(34.75,41.25,N);
% zgrid=linspace(20e3,20e3,N);
% N=50;
% longrid=linspace(139,145,N);
% latgrid=linspace(38,38,N);
% zgrid=linspace(20e3,20e3,N);

% % Large scale griding
% lat=(27:0.05:47);
% lon=(128:0.05:148);
% z=(2:2:100);
% N=max(size(z));
% M=max(size(lon));
% zgrid=repmat(z,M^2,1);
% zgrid=reshape(zgrid,N*M^2,1);
% longrid=repmat(lon,M,1);
% longrid=reshape(longrid,M^2,1);
% longrid=repmat(longrid,N,1);
% latgrid=repmat(lat',M*N,1);

% %Inteprolated
% lonepi=142.369;
% latepi=38.321;
% linelength=7;
% azi=13;
% nsources=30;
% longrid=linspace(lonepi-0.5*linelength*sin(deg2rad(azi)),lonepi+0.5*linelength*sin(deg2rad(azi)),nsources);
% latgrid=linspace(latepi-0.5*linelength*cos(deg2rad(azi)),latepi+0.5*linelength*cos(deg2rad(azi)),nsources);
% zgrid=20e3*ones(size(longrid));


%No. of nodes
N=max(size(longrid));
%N=size(longrid,1);
%Station coordiantes
load(data)
%Select some data
%coseis=tohoku_clean;
% i=find(coseis.process==1);  %Useful data
i=1:1:size(coseis.T,1);
lat=coseis.lat(i)';
lon=coseis.lon(i)';

%Initalize
G=zeros(max(size(i))*3,5,N);
%Loop over inversion nodes and compute GFs
cd(pathgf)
%figure
%hold on
for j=1:N
    display(['Node ' num2str(j) ' of ' num2str(N) ' ...'])
    late = latgrid(j);
    lone = longrid(j);
    depth = zgrid(j);
    %Define current node as origin and compute batch_green('japan','(x,y) coords of stations (in
    %metres)
    %Get station info
    [dist,az] = distance(late,lone,lat,lon);
    dist=deg2km(dist)*1000;
    az=-az+90;
    az=deg2rad(az);
    [e n]=pol2cart(az,dist);
    x=n;
    y=e;
    az=atan2(y,x);
    %scatter(y,x)
    
    
    %MAKE KERNEL MATRIX __________________
    % 1.- G(i,j,k)
    [GZ1 GR1 GT1 GZ2 GR2 GT2 GZ3 GR3 GZ4 GR4 GT4 GZ5 GR5 GT5] =makegreen(fgreen,x',y',depth);
    R=buildrotmat(az,1,1);
    %Get neighboring positions to compute derivatives
    if derivflag==1
        xp=x+dGF;
        xm=x-dGF;
        yp=y+dGF;
        ym=y-dGF;
        %Get new azimuths for rotation from cylindrical to cartesian
        az_xpy=atan2(y,xp);
        az_xpyp=atan2(yp,xp);
        az_xyp =atan2(yp,x);
        az_xmyp=atan2(yp,xm);
        az_xmy=atan2(y,xm);
        az_xmym=atan2(ym,xm);
        az_xym =atan2(ym,x);
        az_xpym=atan2(ym,xp);
        % 2.- G(i+1,j,k)
        [GZ1(:,2) GR1(:,2) GT1(:,2) GZ2(:,2) GR2(:,2) GT2(:,2) GZ3(:,2) GR3(:,2) GZ4(:,2) GR4(:,2) GT4(:,2) ...
            GZ5(:,2) GR5(:,2) GT5(:,2)] =makegreen(fgreen,xp',y',depth);
        R(:,:,2)=buildrotmat(az_xpy,1,1);
        % 3.- G(i-1,j,k)
        [GZ1(:,3) GR1(:,3) GT1(:,3) GZ2(:,3) GR2(:,3) GT2(:,3) GZ3(:,3) GR3(:,3) GZ4(:,3) GR4(:,3) GT4(:,3) ...
            GZ5(:,3) GR5(:,3) GT5(:,3)] =makegreen(fgreen,xm',y',depth);
        R(:,:,3)=buildrotmat(az_xmy,1,1);
        % 4.- G(i,j+1,k)
        [GZ1(:,4) GR1(:,4) GT1(:,4) GZ2(:,4) GR2(:,4) GT2(:,4) GZ3(:,4) GR3(:,4) GZ4(:,4) GR4(:,4) GT4(:,4) ...
            GZ5(:,4) GR5(:,4) GT5(:,4)] =makegreen(fgreen,x',yp',depth);
        R(:,:,4)=buildrotmat(az_xyp,1,1);
        % 5.- G(i,j-1,k)
        [GZ1(:,5) GR1(:,5) GT1(:,5) GZ2(:,5) GR2(:,5) GT2(:,5) GZ3(:,5) GR3(:,5) GZ4(:,5) GR4(:,5) GT4(:,5) ...
            GZ5(:,5) GR5(:,5) GT5(:,5)] =makegreen(fgreen,x',ym',depth);  
        R(:,:,5)=buildrotmat(az_xym,1,1);
        % 6.- G(i,j,k+1)
        [GZ1(:,6) GR1(:,6) GT1(:,6) GZ2(:,6) GR2(:,6) GT2(:,6) GZ3(:,6) GR3(:,6) GZ4(:,6) GR4(:,6) GT4(:,6) ...
            GZ5(:,6) GR5(:,6) GT5(:,6)] =makegreen([fgreen '_z'],x',y',depth);
        R(:,:,6)=buildrotmat(az,1,1);
        % 7.- G(i,j,k+2)
        [GZ1(:,7) GR1(:,7) GT1(:,7) GZ2(:,7) GR2(:,7) GT2(:,7) GZ3(:,7) GR3(:,7) GZ4(:,7) GR4(:,7) GT4(:,7) ...
            GZ5(:,7) GR5(:,7) GT5(:,7)] =makegreen([fgreen '_2z'],x',y',depth);  
        R(:,:,7)=buildrotmat(az,1,1);
        % 8.- G(i+1,j+1,k)
        [GZ1(:,8) GR1(:,8) GT1(:,8) GZ2(:,8) GR2(:,8) GT2(:,8) GZ3(:,8) GR3(:,8) GZ4(:,8) GR4(:,8) GT4(:,8) ...
            GZ5(:,8) GR5(:,8) GT5(:,8)] =makegreen(fgreen,xp',yp',depth);
        R(:,:,8)=buildrotmat(az_xpyp,1,1);
        % 9.- G(i-1,j+1,k)
        [GZ1(:,9) GR1(:,9) GT1(:,9) GZ2(:,9) GR2(:,9) GT2(:,9) GZ3(:,9) GR3(:,9) GZ4(:,9) GR4(:,9) GT4(:,9) ...
            GZ5m_x_yp_z GR5(:,9) GT5(:,9)] =makegreen(fgreen,xm',yp',depth);
        R(:,:,9)=buildrotmat(az_xmyp,1,1);
        % 10.- G(i+1,j-1,k)
        [GZ1(:,10) GR1(:,10) GT1(:,10) GZ2(:,10) GR2(:,10) GT2(:,10) GZ3(:,10) GR3(:,10) GZ4(:,10) GR4(:,10) GT4(:,10) ...
            GZ5(:,10) GR5(:,10) GT5(:,10)] =makegreen(fgreen,xp',ym',depth);
        R(:,:,10)=buildrotmat(az_xpym,1,1);
        % 11.- G(i-1,j-1,k)
        [GZ1(:,11) GR1(:,11) GT1(:,11) GZ2(:,11) GR2(:,11) GT2(:,11) GZ3(:,11) GR3(:,11) GZ4(:,11) GR4(:,11) GT4(:,11) ...
            GZ5(:,11) GR5(:,11) GT5(:,11)] =makegreen(fgreen,xm',ym',depth);
        R(:,:,11)=buildrotmat(az_xmym,1,1);
        % 12.- G(i+1,j,k+1)
        [GZ1(:,12) GR1(:,12) GT1(:,12) GZ2(:,12) GR2(:,12) GT2(:,12) GZ3(:,12) GR3(:,12) GZ4(:,12) GR4(:,12) GT4(:,12) ...
            GZ5(:,12) GR5(:,12) GT5(:,12)] =makegreen([fgreen '_z'],xp',y',depth);
        R(:,:,12)=buildrotmat(az_xpy,1,1);
        % 13.- G(i-1,j,k+1)
        [GZ1(:,13) GR1(:,13) GT1(:,13) GZ2(:,13) GR2(:,13) GT2(:,13) GZ3(:,13) GR3(:,13) GZ4(:,13) GR4(:,13) GT4(:,13) ...
            GZ5(:,13) GR5(:,13) GT5(:,13)] =makegreen([fgreen '_z'],xm',y',depth);
        R(:,:,13)=buildrotmat(az_xmy,1,1);
        % 14.- G(i,j+1,k+1)
        [GZ1(:,14) GR1(:,14) GT1(:,14) GZ2(:,14) GR2(:,14) GT2(:,14) GZ3(:,14) GR3(:,14) GZ4(:,14) GR4(:,14) GT4(:,14) ...
            GZ5(:,14) GR5(:,14) GT5(:,14)] =makegreen([fgreen '_z'],x',yp',depth);
        R(:,:,14)=buildrotmat(az_xyp,1,1);
        % 15.- G(i,j-1,k+1)
        [GZ1(:,15) GR1(:,15) GT1(:,15) GZ2(:,15) GR2(:,15) GT2(:,15) GZ3(:,15) GR3(:,15) GZ4(:,15) GR4(:,15) GT4(:,15) ...
            GZ5(:,15) GR5(:,15) GT5(:,15)] =makegreen([fgreen '_z'],x',ym',depth);
        R(:,:,15)=buildrotmat(az_xym,1,1);
    end
    %Prepare inversion matrices
    display('Preparing matrices for inversion...')
    if derivflag==1
        iter=15;
    else
        iter=1;
    end
    iG=0;
    while iG<iter
        k=1;
        iG=iG+1
        for i=1:size(GR1,1)
            %Vertical component of Gram matrix
            Gt(k,1)=GZ1(i,iG);
            Gt(k,2)=GZ2(i,iG);
            Gt(k,3)=GZ3(i,iG);
            Gt(k,4)=GZ4(i,iG);
            Gt(k,5)=GZ5(i,iG);
            %Radial component of Gram matrix
            Gt(k+1,1)=GR1(i,iG);
            Gt(k+1,2)=GR2(i,iG);
            Gt(k+1,3)=GR3(i,iG);
            Gt(k+1,4)=GR4(i,iG);
            Gt(k+1,5)=GR5(i,iG);
            %Tangential component of Gram matrix
            Gt(k+2,1)=GT1(i,iG);
            Gt(k+2,2)=GT2(i,iG);
            Gt(k+2,3)=0;
            Gt(k+2,4)=GT4(i,iG);
            Gt(k+2,5)=GT5(i,iG);
            k=k+3;
        end
        %Append current matrix
        Gpre(:,:,j,iG)=R(:,:,iG)*Gt;
    end
end
%Now compute the actual derivatives
%G
G=Gpre(:,:,:,1);
if derivflag==1
    %Gx
    G(:,:,:,2)=(Gpre(:,:,:,2)-Gpre(:,:,:,3))/(2*dGF);
    %Gy
    G(:,:,:,3)=(Gpre(:,:,:,4)-Gpre(:,:,:,5))/(2*dGF);
    %Gz/GFs
    G(:,:,:,4)=(-Gpre(:,:,:,7)+4*Gpre(:,:,:,6)-3*Gpre(:,:,:,1))/(2*dGF);
    %Gxx
    G(:,:,:,5)=(Gpre(:,:,:,2)-2*Gpre(:,:,:,1)+Gpre(:,:,:,3))/(dGF^2);
    %Gyy
    G(:,:,:,6)=(Gpre(:,:,:,4)-2*Gpre(:,:,:,1)+Gpre(:,:,:,5))/(dGF^2);
    %Gzz
    G(:,:,:,7)=(Gpre(:,:,:,7)-2*Gpre(:,:,:,6)+Gpre(:,:,:,1))/(dGF^2);
    %Gxy
    G(:,:,:,8)=(Gpre(:,:,:,8)-Gpre(:,:,:,9)-Gpre(:,:,:,10)+Gpre(:,:,:,11))/(4*dGF^2);
    %Gxz
    G(:,:,:,9)=(Gpre(:,:,:,12)-Gpre(:,:,:,2)-Gpre(:,:,:,13)+Gpre(:,:,:,3))/(4*dGF^2);
    %Gyz
    G(:,:,:,10)=(Gpre(:,:,:,14)-Gpre(:,:,:,4)-Gpre(:,:,:,15)+Gpre(:,:,:,5))/(4*dGF^2);
end

cd([pathdata '/GFs'])
grid=vertcat(longrid,vertcat(latgrid,zgrid));
save(fsave,'G','grid')