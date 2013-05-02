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
pathgf='/Users/dmelgarm/Documents/scripts/Fortran/SyntheticGPS/grnfcts'
%Path where data is stored
pathdata='/Users/dmelgarm/Documents/Research/Data/Tohoku'

% Make inversion nodes
cd(pathdata)
[longrid latgrid zgrid]=textread('fault.xyz','%f%f%f');
longrid=longrid';
latgrid=latgrid';
zgrid=zgrid;
zgrid=-zgrid';
% i=find(longrid>134 & longrid<180);
% zgrid=zgrid(i);
% latgrid=latgrid(i);
% longrid=longrid(i);
%i=find(latgrid>40 & latgrid<44);
i=find(zgrid>150);
zgrid=zgrid(i);
latgrid=latgrid(i);
longrid=longrid(i);
display('Creating grid on plate interface...')
% ________________

% % Prism around epicentre
%[longrid latgrid zgrid]=makegrid(-115.28,32.192,10e3,0.1,0.1,2e3,14,14,10);
 %[longrid latgrid zgrid]=makegrid(143.91,41.81,25e3,0.3,0.3,3e3,12,12,10);
% [longrid latgrid zgrid]=makegrid(143.84,42.21,27e3,0.3,0.3,3e3,0,0,0)
% [longrid latgrid zgrid]=makegrid(142.969,38.321,20e3,0.2,0.2,3e3,1,1,1);
% [longrid latgrid zgrid]=makegrid(143.11,37.92,19.5e3,0.2,0.2,3e3,1,1,1);
% i=find(zgrid<150e3);
% longrid=longrid(i);
% latgrid=latgrid(i);
% zgrid=zgrid(i);
% display('Creating prism grid around origin...')
% ____________

% Single sources
% Toki
% longrid=[143.8]; 
% latgrid=[42.1];
% zgrid=[28e3];
% El May
% longrid=[-115.37];
% latgrid=[32.35];
% zgrid=[8e3];
% Tohoku
% longrid=143.11;
% latgrid=37.92;
% zgrid=19.5e3;
% display('No grid created, single point inversion...')
% _____________

%No. of nodes
N=size(longrid,2);
%Station coordiantes
load(data)
lat=coseis.lat;
lon=coseis.lon;
%Initalize
G=zeros(size(coseis.N,1)*3,5,N);
%Loop over inversion nodes and compute GFs
cd(pathgf)
for j=1:N
    display(['Node ' num2str(j) ' of ' num2str(N) ' ...'])
    late = latgrid(j);
    lone = longrid(j);
    depth = zgrid(j);
    %Define current node as origin and compute (x,y) coords of stations (in
    %metres)
    %Get station info
    [dist,az] = distance(late,lone,lat,lon);
    dist=deg2km(dist)*1000;
    az=-az+90;
    az=deg2rad(az);
    [e n]=pol2cart(az,dist);
    x=n;
    y=e;
    
    
    %MAKE KERNEL MATRIX __________________
    [GZ1 GR1 GT1 GZ2 GR2 GT2 GZ3 GR3 GZ4 GR4 GT4 GZ5 GR5 GT5] =makegreen(fgreen,x,y,depth);
    
    %Prepare inversion matrices
    display('Preparing matrices for inversion...')
    k=1;
    for i=1:size(GR1,1)
        %Vertical component of Gram matrix
        Gt(k,1)=GZ1(i);
        Gt(k,2)=GZ2(i);
        Gt(k,3)=GZ3(i);
        Gt(k,4)=GZ4(i);
        Gt(k,5)=GZ5(i);
        %Radial component of Gram matrix
        Gt(k+1,1)=GR1(i);
        Gt(k+1,2)=GR2(i);
        Gt(k+1,3)=GR3(i);
        Gt(k+1,4)=GR4(i);
        Gt(k+1,5)=GR5(i);
        %Tangential component of Gram matrix
        Gt(k+2,1)=GT1(i);
        Gt(k+2,2)=GT2(i);
        Gt(k+2,3)=0;
        Gt(k+2,4)=GT4(i);
        Gt(k+2,5)=GT5(i);
        k=k+3;
    end
    %Append current matrix
    G(:,:,j)=Gt;
end
cd(pathdata)
grid=vertcat(longrid,vertcat(latgrid,zgrid));
save(fsave,'G','grid')