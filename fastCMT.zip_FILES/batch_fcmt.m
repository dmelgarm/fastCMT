function [xyz M L]=batch_fcmt(fsave,Ln)

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
path='/Users/dmelgarm/Documents/Research/Data/Tohoku'
mainpath='/Users/dmelgarm/Documents/scripts/GMT/fastCMT'
data='tohoku.mat';
fgreen='tohokuG.mat';
wflag=1; %Weigh data by preevent noise
wdflag=1; %Weigh flag by distance to centroid (1/d^2)
pflag=0; %Plot data for each epoch and each node
% _________________

%Load data
cd(path)
mkdir(fsave)
load velmod.mat
load tohokuG.mat
load(data)
% _________________

longrid=grid(1,:)';
latgrid=grid(2,:)';
zgrid=grid(3,:)';
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
%Run main loop
for k=1:Niter
    tic
    Gi=G(:,:,k);
    lon=longrid(k)
    lat=latgrid(k)
    z=zgrid(k)
    display([num2str(k) '/' num2str(Niter)])
    [m l d]=fastCMT(coseis,Gi,velmod,[lon,lat z],'psmeca.inp',wflag,pflag,wdflag,Ln);
    cd([path '/' fsave])
    source=[lon lat z];
    sname=['sdisp' num2str(k)];
    save(sname,'source','d');
    M(:,:,:,k)=m;
    L(k,:)=l;
    xyz(k,:)=[lon lat z];
    toc
end
cd([path '/' fsave])
save(fsave,'M','L','xyz')

