function make_psmeca(M,lon,lat,depth)

cd ~/Documents/scripts/GMT/fastCMT
expt=floor(log10(max(max(abs(M)))));
M=M/(1*10^expt);
mrr=M(1,1);
mtt=M(2,2);
mpp=M(3,3);
mrt=M(1,2);
mrp=M(1,3);
mtp=M(2,3);
Mps=[lat lon depth mrr mtt mpp mrt mrp mtp expt];
save('psmeca.inp','Mps','-ascii');
cd ~

