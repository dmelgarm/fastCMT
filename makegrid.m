function [gridlon gridlat gridz]=makegrid(lon,lat,z,dlon,dlat,dz,nlon,nlat,nz)


tlon=lon-floor(nlon/2)*dlon:dlon:lon+floor(nlon/2)*dlon;
tlat=lat-floor(nlat/2)*dlat:dlat:lat+floor(nlat/2)*dlat;
tz=z-floor(nz/2)*dz:dz:z+floor(nz/2)*dz;
i=find(tz>0);
tz=tz(i);
nlat=size(tlat,2);
nlon=size(tlon,2);
nz=size(tz,2);
n=0;
for klon=1:nlon
    klon
    for klat=1:nlat
        for kz=1:nz
            n=n+1;
            gridlon(n)=tlon(klon);
            gridlat(n)=tlat(klat);
            gridz(n)=tz(kz);
        end
    end
end
