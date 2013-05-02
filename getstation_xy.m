function getstation_xy(coseis)

%DMM 10/2011

%Turn lat lon station coordinates into (x,y)

late=42.21;
lone=143.84;
i=1:1:size(coseis.T,1);
lat=coseis.lat(i);
lon=coseis.lon(i);
[dist,az] = distance(late,lone,lat,lon);
dist=deg2km(dist)*1000;
az=-az+90;
az=deg2rad(az);
[e n]=pol2cart(az,dist);
x=n;
y=e;
cd('/diego-local/scripts/Fortran/ED')
delete station.xy
fid = fopen('station.xy', 'w');
for k=1:length(x)
    sta=['( ' num2str(x(k)) ', ' num2str(y(k)) '),'];
    fprintf(fid, '( %4.3e, %4.3e),\n', x(k), y(k));
    clear sta
end
close(fid)
%Write to file
cd('/diego-local/scripts/Fortran/ED')
delete station.xy
fid = fopen('station.xy', 'w');
fprintf(fid, '%s\n', sta);
close(fid)