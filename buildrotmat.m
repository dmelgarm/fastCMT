function R=buildrotmat(azi,direction);

%DMM 01/2011
%
%EDGRN code outptus Green functions in a local cilindircal frame, here we
%build a rotation matrix from local cilindircal to local cartesian given a
%source-station azimuth (from north).
%The matrix is built in blocs of 3x3 where each block has it's own azimuth
%corresponding to one station-event pair.

%direction -  1 is convert from z,r,theta to z,x,y
%             2 is convert from z,x,y to z,r,theta


n=max(size(azi,1),size(azi,2));
R=zeros(n*3,n*3);
if direction==1
    for k=1:n
        r=[1 0 0;0 cos(azi(k)) -sin(azi(k));0 sin(azi(k)) cos(azi(k))];
        j=3*(k-1)+1;
        R(j:j+2,j:j+2)=r;
    end
else
    for k=1:n
        r=[1 0 0;0 cos(azi(k)) sin(azi(k));0 -sin(azi(k)) cos(azi(k))];
        j=3*(k-1)+1;
        R(j:j+2,j:j+2)=r;
    end
end