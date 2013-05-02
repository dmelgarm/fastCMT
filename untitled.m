function [longrid latgrid zgrid]=getgrid(lon,lat,z,azi,length,nsource)

%DMM 03/2012


cd /diego-local/Research/Data/Tohoku/DenseGF
[longrid latgrid zgrid]=textread('japan_aria_grid','%f%f%f');
%Get grid params
dlat
