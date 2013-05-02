function [xyz M L]=batch_fcmt2(fsave,Ln)

%DMM 04/2011
%
% Perform CMT inversion at multiple points using the fastCMT routine
%
% IN
% fsave - file name to save data
% Ln - Norm to be used in the inversions, 1 for L1, 2 for L2
%         
% OUT
% xyz - Contains coordinates of grid points
% M - 4 dimensional variable containing momoent tensors of all grid points at all time steps
% L - Contains misfits (variance reduction) for each grid point at each time step

warning off
format compact
grid=[];

%Run parameters
%path='/diego-local/Research/Data/Tohoku'
%path='/diego-local/Research/Data/Escenario2011'
path='/diego-local/Research/Data/Toki'
%path='/diego-local/Research/Data/El Mayor'
mainpath='/diego-local/scripts/GMT/fastCMT'
%linelength=7; %In degrees
%nsources=30;
%  linelength=5; %In degrees
%  nsources=21;
 linelength=0; %In degrees
 nsources=1;

% %Epicentre coords (This is used in SOPAC gridsearch Tohoku)
% lonepi=141.8:0.25:143.8;
% latepi=37:0.25:39.0;
% depthepi=10e3:10e3:40e3;
% %Line griding
% azi=15:15:180;

% %Epicentre coords (This is used in SINGLE CMT Tohoku)
% lonepi=140.4:0.2:144.4;
% latepi=36.2:0.2:40.4;
% depthepi=10e3:5e3:60e3;
% azi=0;

% %Epicentre coords (This is used in SOPAC gridsearch Toki)
% lonepi=143.2:0.2:144.8;
% latepi=41.4:0.2:43;cd /diego-local/Research/Data/El' Mayor'/
% depthepi=20e3:10e3:50e3;
% %Line griding
% azi=15:15:180;

%Epicentre coords (This is used in SINGLE CMT Toki)
lonepi=142.0:0.2:146;
latepi=40:0.2:44;
depthepi=10e3:5e3:60e3;
azi=0;

% lonepi=142;  %(This is used in SOPAC individual Tohoku)
% latepi=37.2;
% depthepi=30e3;
% azi=15;
%cd /diego-local/Research/Data/El' Mayor'/
% %Epicentre coords (This is used in SOPAC Inidividual [longrid latgrid zgrid inode]=getgrid(lon,lat,z,azi,linelength,nsources,grid);Toki)
% lonepi=144;
% latepi=42.2;
% depthepi=30e3;
% Line griding
% azi=49;

% Epicentre coords (This is used in SOPAC Inidividual El Mayor)
% lonepi=-115.68;
% latepi=32.39;
% depthepi=10e3;
% azi=138;


% %Epicentre coords (This is used in SOPAC gridsearch El Mayor)
% lonepi=-116.4:0.2:-114.8;
% latepi=31.4:0.2:33;
% depthepi=5e3:5e3:25e3;
% %Line griding
% azi=15:15:180;


% %Epicentre coords  MX2011
% lonepi=-101:0.25:-99;
% latepi=17:0.25:18;
% depthepi=5e3:10e3:30e3;
% %Line griding
% azi=105;

% %Epicentre coords  MX2011
% lonepi=-99.5
% latepi=16.75
% depthepi=10e3;
% %Line griding
% azi=105;


% %Epicentre coords
% lonepi=141.3:0.2:142.3;
% latepi=37:0.2:38.0;
% depthepi=5e3:2e3:15e3;
% %Line griding
% azi=30:2:60;

% %Epicentre coords
% lonepi=141.869:0.2:142.869;%142.369;
% latepi=37.821:0.2:38.821;%38.321;
% depthepi=10e3:5e3:60e3;
% lonepi=142.369;
% latepi=38.321;
% depthepi=10e3:5e3:60e3;
% %Line griding
% azi=[13];

% %Epicentre coords
% lonepi=142.5;
% latepi=38.2;
% depthepi=20e3;
% %Line griding
% %azi=5:5:180;
% azi=15;

%Compute GFs on the fly
% lonepi=142.369;
% latepi=38.321;
% depthepi=20e3;
% %Line griding
% azi=[13];
% %MNodes
% longrid=linspace(lonepi-0.5*linelength*sin(deg2rad(azi)),lonepi+0.5*linelength*sin(deg2rad(azi)),nsources);
% latgrid=linspace(latepi-0.5*linelength*cos(deg2rad(azi)),latepi+0.5*linelength*cos(deg2rad(azi)),nsources);
% % a=0;



%Make epi and azi vector
i=0;
for k1=1:max(size(lonepi))
    for k2=1:max(size(latepi))
        for k3=1:max(size(depthepi))
            for k4=1:max(size(azi))
                i=i+1;
                lonall(i)=lonepi(k1);
                latall(i)=latepi(k2);
                depthall(i)=depthepi(k3);
                aziall(i)=azi(k4);
            end
        end
    end
end
zgrid=0; %Otherwise MATLAB wants to use it as a function
% row=1055;
% aziall=aziall(row);
% depthall=depthall(row);
% latall=latall(row);
% lonall=lonall(row);
nlines=max(size(aziall));

%data='toki300.mat';
data='rapid_coseis.mat';
%data='rapid_coseis_allgps.mat';
%data='tohoku_ARIA.mat'; %ARIA offsets
%data='tohoku_gps_all.mat';
%data='tohoku_multi.mat';
%data='tohoku_finite.mat';


%fgreen='50sta_fgreen.mat';
%fgreen='slab1.0_GF.mat';
%fgreen='toki_prism_small_GF.mat';
%fgreen='mx_fplaneGF.mat';
%fgreen='gcmt_GF.mat';
%fgreen='GFmulti.mat'
%fgreen='GFderiv_ARIA.mat'
%fgreen='GFprism_ARIA.mat'
%fgreen='GFmultifine_ARIA.mat'
%fgreen='tohoku_USGS_epi.mat'
%fgreen='japan_ariaGF_'
%fgreen='japan_sopacGF_'
fgreen='japan_northGF_'
%fgreen='elmayGF_'
%fgreen='mx11_GFGF_'
%fgreen='GF_aria_finite';
%fgreen='rapid_kalmanGF_';

wflag=1; %Weigh data by preevent noise
wdflag=1; %Weigh flag by distance to centroid (1/d^2)
dcflag=5; %Force double couple? 5
pflag=0; %Plot data for each epoch and each node
multiflag=0;  %Multiple CMT inversion?
secondflag=0; %2nd order MT inversion?
regflag=0;  %Regularize?
tikh=2;
nlambda=100;
lambda=logspace(-20,50,nlambda);

% _________________

%Load data
display('Loading data and Green functions...')
cd(path)
mkdir(fsave)
load velmod.mat
%load velmod_4layer.mat
load(data)
nsta=1:1:length(coseis.N);


% _________________

% longrid=grid(1,:)';
% latgrid=grid(2,:)';
% zgrid=grid(3,:)';
% longrid=143.11;
% latgrid=37.92;
% zgrid=19500;
% igrd(1)=find(latgrid==36.00 & longrid==142.00 & zgrid==20e3);
% igrd(2)=find(latgrid==36.75 & longrid==142.25 & zgrid==20e3);
% igrd(3)=find(latgrid==37.25 & longrid==142.compact flash50 & zgrid==20e3);
% igrd(4)=find(latgrid==38.00 & longrid==142.75 & zgrid==20e3);
% igrd(5)=find(latgrid==38.75 & longrid==143.00 & zgrid==20e3);
% igrd(6)=find(latgrid==39.50 & longrid==143cd /diego-local/Research/Data/El' Mayor'/.25 & zgrid==20e3);
% igrd(7)=find(latgrid==40.00 & longrid==143.50 & zgrid==20e3);
% igrd=find(latgrid==38 & longrid==143 & zgrid==20e3);
% latgrid=latgrid(igrd);
% longrid=longrid(igrd);
% zgrid=zgrid(igrd);
% G=G(:,:,igrd);



%Get solutions
k=0;
xyz=[0 0 0];
maxdata=0;
display(['L' num2str(Ln) ' norm inversion selected...'])
%Which data are useful
%i=find(tohoku.process==1);


%Extract only relevant rows of Gram matrix
% ii1=nsta*3-2;
% ii2=nsta*3-1;
% ii3=nsta*3;
% j1=interleave(ii1,ii2);92
% j2=interleave(ii2,ii3);
% j=unique(interleave(j1,j2));
% G=G(j,:,:,:);
%Remove Stations that are NaN's
inan=find(~isnan(coseis.N));
%coseis.gname=coseis.gname(inan);
coseis.N=coseis.N(inan,1);
coseis.E=coseis.E(inan,1);
coseis.U=coseis.U(inan,1);
coseis.T=coseis.T(inan,1);
% coseis.Ngps=coseis.Ngps(inan,:);
% coseis.Egps=coseis.Egps(inan,:);
% coseis.Ugps=coseis.Ugps(inan,:);
coseis.stdn=coseis.stdn(inan);
coseis.stde=coseis.stde(inan);
coseis.stdu=coseis.stdu(inan);
coseis.lat=coseis.lat(inan);
coseis.lon=coseis.lon(inan);
inan1=inan*3-2;
inan2=inan*3-1;
inan3=inan*3;
jnan1=interleave(inan1,inan2);
jnan2=interleave(inan2,inan3);
jnan=unique(interleave(jnan1,jnan2));
% G=G(jnan,:,:,:);
if multiflag==1
    Niter=nlines;
else
    Niter=length(lonall);
end

%Run main loop
for k=1:Niter
    tic
    display([num2str(k) '/' num2str(Niter)])
    if multiflag==0
        cd([path '/denseGF'])
        load nodes_grid.mat
        lon=lonall(k);
        lat=latall(k);
        z=depthall(k);
        [longrid latgrid zgrid inode]=getgrid(lon,lat,z,azi,linelength,nsources,grid);
        Gtemp=load([fgreen num2str(inode)]);
        Gi=Gtemp(jnan,:);
        [m moment l d]=fastCMT(coseis,Gi,velmod,[longrid,latgrid zgrid],'psmeca.inp',wflag,pflag,wdflag,Ln,dcflag,secondflag);
    else
        %Load GFs
        cd([path '/denseGF'])
        load nodes_grid.mat
        [longrid latgrid zgrid inode]=getgrid(lonall(k),latall(k),depthall(k),aziall(k),linelength,nsources,grid);
%         cd denseGF
%         load rapid_kalman_grid.mat
%         longrid=grid(1,:);
%         latgrid=grid(2,:);
%         zgrid=grid(3,:);
%         inode=1:length(longrid);
        nnode=max(size(inode))
        for j=1:nnode
            Gtemp=load([fgreen num2str(inode(j))]);
            G(:,:,j)=Gtemp(jnan,:);
        end
%         %Use interpolated GFs
%         cd('/diego-local/Research/Data/Tohoku/GFs')
%         load GF_sopac_interp.mat
%         G=G(jnan,:,:);
%         longrid=grid(1,:)';
%         latgrid=grid(2,:)';
%         zrid=grid(3,:)';
%         %
        Gi=G;
        clear G
        lon=longrid;
        lat=latgrid;
        z=zgrid;
        aziall(k)
        if (max(size(coseis.N))*3)~=max(size(Gi))
            display('FATAL ERROR: Number of stations in data is NOT equal to number of stations in Green functions...')
            display('Exiting...')
            return
        end
        cd(path)
        [m mall moment l d lambda_corner(k) Lcurve(:,1) Lcurve(:,2) Lcurve(:,3)]=mfastCMT(coseis,Gi,velmod,[lon,lat z],...
            'psmeca.inp',wflag,pflag,wdflag,Ln,dcflag,regflag,lambda,tikh);
        M0=norm(m,'fro')/sqrt(2); %Scalar moment
    end
    cd([path '/' fsave])
    source=[lon lat z];
    sname=['sdisp' num2str(k)];
    sourcename=['source' num2str(k)];
    xyzname=['xyz' num2str(k)];
    momentname=['moment' num2str(k)];
    mname=['mall' num2str(k)];
    lcurvename=['Lcurve' num2str(k)];
    M(:,:,:,k)=m;
    L(k,:)=l;
    if multiflag==1
        xyz=[lon lat z];
        mavg=m;
        mall=m;
        save(lcurvename,'Lcurve');
        save(mname,'mall');
    else
        xyz=[longrid' latgrid' zgrid'];
    end
    save(sname,'d');
    save(sourcename,'source');
    save(xyzname,'xyz');
    save(momentname,'moment');
    toc
end
toc
cd([path '/' fsave])
save(fsave,'M','L','xyz')
if multiflag==1
    save('lambda','lambda_corner')
end
display([num2str(length(nsta)) ' stations used in ivnersion.'])
grd=[longrid latgrid zgrid];
cd(mainpath)
save('grdsearch.xyz','grd','-ascii')
cd(path)