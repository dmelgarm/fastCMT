function xy2latlon(x,y,lat,lon)
    
    az=90-rad2deg(angle(x+i*y))
    dist=km2deg(sqrt(x.^2+y.^2));
    [latout,lonout] = reckon(lat,lon,dist,az)
    scatter(lonout,latout),hold on,scatter(lon,lat),grid on,axis equal
    coords=[lonout' latout'];
    coords=[coords;coords(1,:)];
    cd /diego-local/scripts/GMT/fastCMT/
    fid=fopen('rectangular_fault','w');
    fprintf(fid, '%8.4f\t%8.4f\n', coords');
    fclose(fid);
    
