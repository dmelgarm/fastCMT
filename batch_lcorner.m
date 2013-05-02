function batch_lcorner

cd ~/Research/Data/Tohoku/fastCMTtohoku_sopac_gridsearch_1deriv
load fastCMTtohoku_sopac_gridsearch_1deriv
% cd /diego-local/Research/Data/Escenario2011/fastCMTmx11_gridsearch
% load fastCMTmx11_gridsearch
% cd /diego-local/Research/Data/Toki/fastCMTtoki_gridsearch
% load fastCMTtoki_gridsearch
% cd /diego-local/Research/Data/El' Mayor'/fastCMTelmay_gridsearch
% load fastCMTelmay_gridsearch

ha = tight_subplot(2, 1, 0.1, 0.2, 0.37)
N=max(size(L));
axes(ha(1))
for k=1:N
    k
    load(['Lcurve' num2str(k) '.mat'])
    [ppx,kappa,lambda(k),misfit(k),msize(k)] = l_corner(Lcurve(:,2),Lcurve(:,1),Lcurve(:,3));
    loglog(Lcurve(:,2),Lcurve(:,1),'b');
    grid on
    hold on
    scatter(misfit(k),msize(k),'x')
end
xlabel('||WGm-Wd||_1')
ylabel('||Lm||_1')
xi=log10(logspace(log10(min(lambda)),log10(max(lambda)),1000));
[f] = ksdensity(log10(lambda),xi);
axes(ha(2))
plot(xi,f)
grid on
xlabel('log(\lambda)')
ylabel('PDF')
