function [loninterp latinterp]=bestline(longrid,latgrid)

%DMM 04/2012
%
%Get best fitting line for lat lon coordinates sampled at number of point
%sources

y=latgrid;
A=[longrid ones(size(longrid))];
x=lsqlin(A,y);
minterp=x(1);
binterp=x(2);
N=max(size(longrid));
loninterp=linspace(longrid(1),longrid(end),N);
latinterp=loninterp*minterp+binterp;