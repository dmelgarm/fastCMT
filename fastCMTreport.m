function fastCMTreport(m,lonmeca,latmeca,zmeca,L,nsta)

%DMM 10/2011
% 
% Make CMT final solution report files for GMT plotting
% 
% IN:
% M ~ Array of momeont tensors
% 


%Basic stuff
time=180;
var=max(L);
zmeca=zmeca/1000;

%get eigen values, vectors, plunges and azimuths
[teig peig beig tazim tplg pazim pplg bazim bplg]=mom2eigs(m);
%Get exponent and full MT components
mt=mcart2mrp(m);
scale=fix(max(max(log10(abs(mt)))));
scale=1*10^scale;
mt=mt/scale;
mrr=mt(1,1);
mtt=mt(2,2);
mpp=mt(3,3);
mrt=mt(1,2);
mrp=mt(1,3);
mtp=mt(2,3);
teig=teig/scale;
peig=peig/scale;
beig=beig/scale;
scale=fix(log10(scale));
%Make double couple get nodal planes and moment
[mdc e]=mom2dc(m);
Mo=norm(mdc,'fro')/sqrt(2);
Mw=0.67*(log10(Mo)-9.1);
[NP1 NP2 res]=mom2sdr(mdc);


%Write files

make_psmeca_report(m*1e7,lonmeca,latmeca,zmeca*1000);
cd('~/scripts/GMT/fastCMT')
delete fatsCMT.tbl
fid = fopen('fastCMT.tbl', 'w');
fprintf(fid, '%5.2f\n', latmeca);
fprintf(fid, '%5.2f\n', lonmeca);
fprintf(fid, '%3.0f\n', zmeca);
fprintf(fid, '%3.0f\n', time);
fprintf(fid, '%2.1f\n', Mw);
fprintf(fid, '%2.0f\n', scale);
fprintf(fid, '%4.3f\n', mrr);
fprintf(fid, '%4.3f\n', mtt);
fprintf(fid, '%4.3f\n', mpp);
fprintf(fid, '%4.3f\n', mrt);
fprintf(fid, '%4.3f\n', mrp);
fprintf(fid, '%4.3f\n', mtp);
fprintf(fid, '%4.3f\n', teig);
fprintf(fid, '%3.0f\n', tazim);
fprintf(fid, '%2.0f\n', tplg);
fprintf(fid, '%4.3f\n', peig);
fprintf(fid, '%3.0f\n', pazim);
fprintf(fid, '%2.0f\n', pplg);
fprintf(fid, '%4.3f\n', beig);
fprintf(fid, '%3.0f\n', bazim);
fprintf(fid, '%2.0f\n', bplg);
fprintf(fid, '%4.3e\n', Mo);
fprintf(fid, '%3.0f\n', NP1(1));
fprintf(fid, '%3.0f\n', NP1(2));
fprintf(fid, '%3.0f\n', NP1(3));
fprintf(fid, '%3.0f\n', NP2(1));
fprintf(fid, '%3.0f\n', NP2(2));
fprintf(fid, '%3.0f\n', NP2(3));
fprintf(fid, '%4.0f\n', nsta);
fprintf(fid, '%4.2f\n', var);
fclose(fid)


%save('fastCMT.tbl','table','-ascii');
