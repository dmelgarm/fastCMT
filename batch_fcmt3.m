function [xyz M L]=batch_fcmt3(fsave,Ln)

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

%Run parameters
path='/diego-local/Research/Data/El Mayor'
mainpath='/diego-local/scripts/GMT/fastCMT'
%nsta=get_max_coseis;
%nsta=1:1:355;
nsta=1:1:50;
data='d_all_reduced_300.mat';
%data='synthFinite.mat';
fgreen='50sta_fgreen.mat';
wflag=1; %Weigh data by preevent noise
wdflag=1; %Weigh flag by distance to centroid (1/d^2)
dcflag=5; %Force double couple? 5 for no 4 for yes
pflag=0; %Plot data for each epoch and each node
decimateflag=0;  %Downsample data

% _________________

%Load data
display('Loading data and Green functions...')
cd(path)
mkdir(fsave)
load velmod.mat
%load velmod_4layer.mat
load(fgreen)
load(data)
% _________________

longrid=grid(1,:)';
latgrid=grid(2,:)';
zgrid=grid(3,:)';
% longrid=143.11;
% latgrid=37.92;
% zgrid=19500;
grd=[longrid latgrid zgrid];
cd(mainpath)
save('grdsearch.xyz','grd','-ascii')
cd(path)

%Structural prelims
Niter=max(size(latgrid));
%Get solutions
k=0;
xyz=[0 0 0];
maxdata=0;
display(['L' num2str(Ln) ' norm inversion selected...'])
%Which data are useful
%i=find(tohoku.process==1);
if decimateflag==1
    display('Decimating data...')
    R=100;
    for k=1:max(size(tohoku_clean.lat))
        coseis.N(:,k)=downsample(tohoku_clean.nmagps(:,k),R);
        coseis.E(:,k)=downsample(tohoku_clean.emagps(:,k),R);
        coseis.U(:,k)=downsample(tohoku_clean.umagps(:,k),R);
        coseis.T(:,k)=rnd(downsample(tohoku_clean.tmagps(:,k),R));
    end
    coseis.lat=tohoku_clean.lat;
    coseis.lon=tohoku_clean.lon;
    coseis.stdn=tohoku_clean.stdngps;
    coseis.stde=tohoku_clean.stdegps;
    coseis.stdu=tohoku_clean.stdugps;
    clear tohoku_clean
else
%     coseis.N=tohoku_clean.nmagps;
%     coseis.E=tohoku_clean.emagps;
%     coseis.U=tohoku_clean.umagps;
%     coseis.T=tohoku_clean.tmagps;
%     coseis.lat=tohoku_clean.lat;
%     coseis.lon=tohoku_clean.lon;
%     coseis.stdn=tohoku_clean.stdngps;
%     coseis.stde=tohoku_clean.stdegps;
%     coseis.stdu=tohoku_clean.stdugps;
%     clear tohoku_clean
end

%Extract only relevant rows of Gram matrix
ii1=nsta*3-2;
ii2=nsta*3-1;
ii3=nsta*3;
j1=interleave(ii1,ii2);
j2=interleave(ii2,ii3);
j=unique(interleave(j1,j2));
G=G(j,:,:);


%Run main loop
for k=1:Niter
    if k==NaN
        pflag=1
    else
        pflag=0
    end
    tic
    Gi=G(:,:,k);
    lon=longrid(k)
    lat=latgrid(k)
    z=zgrid(k)
    display([num2str(k) '/' num2str(Niter)])
    [m l d]=fastCMT3(coseis,Gi,velmod,[lon,lat z],'psmeca.inp',wflag,pflag,wdflag,Ln,dcflag);
    cd([path '/' fsave])
    source=[lon lat z];
    sname=['sdisp' num2str(k)];
    save(sname,'source','d');
    M(:,:,:,k)=m;
    L(k,:)=l;
    xyz(k,:)=[lon lat z];
    toc
    if k==NaN
        pause
    end
end
cd([path '/' fsave])
save(fsave,'M','L','xyz')
display([num2str(length(nsta)) ' stations used in ivnersion.'])
