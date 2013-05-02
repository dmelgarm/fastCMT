function plot_moments(mkal,mgps,lat)

for k=1:length(mgps)
    Mgps(k)=norm(mgps(:,:,k),'fro')/sqrt(2);
    Mkal(k)=norm(mkal(:,:,k),'fro')/sqrt(2);
end

ha = tight_subplot(1, 1, 0, 0.3, 0.3)
axes(ha(1))
plot(lat,Mkal/10^21,lat,Mgps/10^21)
xlabel('North Latitude (°)')
ylabel('Moment (Nmx10^2^1)')
legend('Kalman','GPS')