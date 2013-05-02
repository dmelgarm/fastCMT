function update_toki

cd /diego-local/Research/Data/Toki
load toki_kal.mat
coseisk=coseis;
load tokiall.mat
i=coseisk.gps;
for k=1:length(i);
    k
    if k==4
        a=0;
    end
    sta=i(k);
    j=find(coseis.gpssta==sta);
    coseisk.N(k,:)=coseis.N(j,:);
    coseisk.E(k,:)=coseis.E(j,:);
    coseisk.U(k,:)=coseis.U(j,:);
end
clear coseis
coseis=coseisk;
save('toki_kal.mat','coseis');