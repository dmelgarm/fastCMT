function [GZ1 GR1 GT1 GZ2 GR2 GT2 GZ3 GR3 GZ4 GR4 GT4 GZ5 GR5 GT5] = makegreen(f,XS,YS,ze)

%DMM 01/2011
%
%Read GFs obtained through EGRNF and output GF values at station
%coordinates. EDGRN produces GFs at discrete radial distances from the 
%source, here we use spline interpolation to obtain the value of the GF at
%a given point and then weigh by a trig. function of the azimuth. This last
%part is critical. For a source-stationn azimuth angle 'azi' from north one
%must compute the radial(r), transverse(t) and vertical (z) component of
%the Green function. Each station event pair has 5 GFs, ecah GF is mono
%valued and must be weighed as follows:
%
%G1=(G1z,G1r,G1t) is pure strike slip
%   G1z=Gz*sin(2*azi)
%   G1r=Gr*sin(2*azi)
%   G1t=Gt*cos(2*azi)
%
%G2=(G2z,G2r,G2t) is pure dip slip
%   G2z=Gz*cos(azi)
%   G2r=Gr*cos(azi)
%   G2t=Gt*sin(azi)
%
%G3=(G3z,G3r,0) is compensated linear vector dipole (CLVD) and is axis-symmetric
%   G3z=Gz
%   G3r=Gr
%   G3t=0
%
%G4=(G4z,G4r,G4t) is N45E strike slip
%   G4z=Gz*cos(2*azi)
%   G4r=Gr*cos(2*azi)
%   G4t=Gt*-sin(2*azi)
%
%G1=(G5z,G5r,G5t) is 45 degree dip dip-slip
%   G5z=Gz*sin(azi)
%   G5r=Gr*sin(azi)
%   G5t=Gt*-cos(azi)



warning('off','all')
display('Interpolating GFs...')


%Read file from EDGRN
fid=fopen([f '.ss']);
%first 12 lines are headers, 13th line contains info on file structure,
%keep it!
for k=1:12
    fstruc=fgets(fid);
end
fstruc=str2num(fstruc);

%Get useful parameters
nr=fstruc(1);   %No. of radial distances
d1=fstruc(2);   %Start and stop radii
d2=fstruc(3);
nzs=fstruc(4); %No. of depth slices
z1=fstruc(5);  %Start and stop depths
z2=fstruc(6);
z0=fstruc(7);  %Observation depth
fclose(fid);
dz=(z2-z1)/(nzs-1);  %Depth discretization
dr=(d2-d1)/(nr-1);  %Range disretization

%Get station azimuths and distances to origin
d=sqrt(XS.^2+YS.^2);
i=find(d<d2);  %Keep only stations in range
ista=i;  %Index of stations I'm keeping
d=d(i);
XS=XS(i);
YS=YS(i);
az=atan2(YS,XS);  %get azimuth
nsta=max([size(XS,1) size(XS,2)]);   %stations left
AZ=repmat(az,1,nr);   %For multiplying times GF matrix

%Get green functions
[Gsstemp Gdstemp Gcltemp]=rgreen(f);
%Keep dispalcements only
Gsstemp=Gsstemp(:,1:3);
Gdstemp=Gdstemp(:,1:3);
Gcltemp=Gcltemp(:,1:2);
%Re-tile into 3D array
iread=1:1:nr;
Gss=zeros(nr,3,nzs);
Gds=zeros(nr,3,nzs);
Gcl=zeros(nr,2,nzs);
for k=1:nzs
    i=iread+(k-1)*nr;
    Gss(:,:,k)=Gsstemp(i,:);
    Gds(:,:,k)=Gdstemp(i,:);
    Gcl(:,:,k)=Gcltemp(i,:);
end
%Determine position in array for interpolation of GFs for each station
iz=floor((ze-z1)/dz)+1;  %How many discretization points from obs depth to source depth
dzs=(ze-(z1+dz*(iz-1)))/dz;   %Discretized source depth
idis=floor((d-d1)/dr)+1;
ddis=(d-(d1+dr*(idis-1)))/dr;
%Interpolation cosntants
w00=(1-ddis).*(1-dzs);
w10=ddis.*(1-dzs);
w01=(1-ddis).*dzs;
w11=ddis.*dzs;
%Interpolate GFs at station coordinates
%z component
Gssz=w00.*Gss(idis,1,iz)+w10.*Gss(idis+1,1,iz)+w01.*Gss(idis,1,iz+1)+...
    w11.*Gss(idis+1,1,iz+1);
Gdsz=w00.*Gds(idis,1,iz)+w10.*Gds(idis+1,1,iz)+w01.*Gds(idis,1,iz+1)+...
    w11.*Gds(idis+1,1,iz+1);
Gclz=w00.*Gcl(idis,1,iz)+w10.*Gcl(idis+1,1,iz)+w01.*Gcl(idis,1,iz+1)+...
    w11.*Gcl(idis+1,1,iz+1);
%r component
Gssr=w00.*Gss(idis,2,iz)+w10.*Gss(idis+1,2,iz)+w01.*Gss(idis,2,iz+1)+...
    w11.*Gss(idis+1,2,iz+1);
Gdsr=w00.*Gds(idis,2,iz)+w10.*Gds(idis+1,2,iz)+w01.*Gds(idis,2,iz+1)+...
    w11.*Gds(idis+1,2,iz+1);
Gclr=w00.*Gcl(idis,2,iz)+w10.*Gcl(idis+1,2,iz)+w01.*Gcl(idis,2,iz+1)+...
    w11.*Gcl(idis+1,2,iz+1);
%t component
Gsst=w00.*Gss(idis,3,iz)+w10.*Gss(idis+1,3,iz)+w01.*Gss(idis,3,iz+1)+...
    w11.*Gss(idis+1,3,iz+1);
Gdst=w00.*Gds(idis,3,iz)+w10.*Gds(idis+1,3,iz)+w01.*Gds(idis,3,iz+1)+...
    w11.*Gds(idis+1,3,iz+1);
%Apply azimuthal scaling to GFs
%z comp.
GZ1=Gssz.*sin(2*az);
GZ2=Gdsz.*cos(az);
GZ3=Gclz;
GZ4=Gssz.*cos(2*az);
GZ5=Gdsz.*sin(az);
%r comp.
GR1=Gssr.*sin(2*az);
GR2=Gdsr.*cos(az);
GR3=Gclr;
GR4=Gssr.*cos(2*az);
GR5=Gdsr.*sin(az);
%z comp.
GT1=Gsst.*cos(2*az);
GT2=Gdst.*sin(az);
GT4=Gsst.*(-sin(2*az));
GT5=Gdst.*(-cos(az));


display(['GFs interpolated up to threshold distance ' num2str(d2/1000) 'km'])



