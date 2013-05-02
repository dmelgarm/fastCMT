function plot_edcmp(f,flag)
%DMM 10/2011
%
% flag = 0 for point source
% flag = 1 for finite fault


display('Reading displacements from EDCMP...')
cd '/diego-local/scripts/Fortran/ED/outdata'
fid=fopen([f '.disp']);
%Ignore header info (3 lines)
for k=1:3
    t=fgets(fid);
end
D=fscanf(fid,'%f%f%f%f%f',[5 inf]);
fclose(fid);
D=D';
y=D(:,1);
x=D(:,2);
uy=D(:,3);
ux=D(:,4);
uz=D(:,5);
clear D
%get lat's and lon's
cd /diego-local/Research/Data/Toki
load toki300.mat
figure
quiver(coseis.lon,coseis.lat,ux,uy,5)
grid on
axis equal
%Write to file
coseis.N=uy;
coseis.E=ux;
coseis.U=uz;
coseis.T=ones(size(coseis.N));
coseis.stdn=ones(size(coseis.N));
coseis.stde=ones(size(coseis.N));
coseis.stdu=ones(size(coseis.N));
if flag==0;
    save('synthCMT.mat','coseis')
else
    save('synthFinite.mat','coseis')
end

