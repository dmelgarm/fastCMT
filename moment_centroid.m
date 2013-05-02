function [lonmeca latmeca]=moment_centroid(moment)

N=100;
%moment=flipud(moment);
phi=atan((moment(end,2)-moment(1,2))/(moment(end,1)-moment(1,1)));
dist=sqrt((moment(:,1)-moment(1,1)).^2+(moment(:,2)-moment(1,2)).^2);
d=linspace(dist(1),dist(end),N);
m=interp1(dist,moment(:,3),d);
M=trapz(d,m);
k=2;
while k<N
    Mtest=trapz(d(1:k),m(1:k));
    if Mtest>M/2 %centroid found
        dcentroid=0.5*(d(k-1)+d(k));
        k=N+1;
    end
    k=k+1;
end
lonmeca=moment(1,1)+dcentroid*cos(phi);
latmeca=moment(1,2)+dcentroid*sin(phi);