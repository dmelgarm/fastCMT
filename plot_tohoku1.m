function plot_tohoku1(tohoku,sta)

%DMM 10/2011
%
%Plot one station


i=find(tohoku.gname==sta);
figure
an=tohoku.an(:,i)/1000;
ae=tohoku.ae(:,i)/1000;
au=tohoku.au(:,i)/1000;
ta=tohoku.ta(:,i);
dkn=tohoku.dkn(:,i);
dke=tohoku.dke(:,i);
dku=tohoku.dku(:,i);
vkn=tohoku.vkn(:,i);
vke=tohoku.vke(:,i);
vku=tohoku.vku(:,i);
tk=tohoku.tk(:,i);
ngps=tohoku.ngps(:,i);
egps=tohoku.egps(:,i);
ugps=tohoku.ugps(:,i);
tgps=tohoku.tgps(:,i);

subplot(3,3,1)
plot(ta,an,'LineWidth',2)
grid on
ylabel('Accel. (m/s^2)','FontSize',18)
title('North','FontSize',18)
xlim([5 250])
set(gca,'FontSize',16)
subplot(3,3,2)
plot(ta,ae,'LineWidth',2)
grid on
title('East','FontSize',18)
xlim([5 250])
set(gca,'FontSize',16)
subplot(3,3,3)
plot(ta,au,'LineWidth',2)
grid on
title('Up','FontSize',18)
xlim([5 250])
set(gca,'FontSize',16)

subplot(3,3,4)
plot(tk,vkn,'LineWidth',2)
grid on
ylabel('Vel. (m/s)','FontSize',18)
xlim([5 250])
set(gca,'FontSize',16)
subplot(3,3,5)
plot(tk,vke,'LineWidth',2)
grid on
xlim([5 250])
set(gca,'FontSize',16)
subplot(3,3,6)
plot(tk,vku,'LineWidth',2)
grid on
xlim([5 250])
set(gca,'FontSize',16)

subplot(3,3,7)
plot(tgps,ngps,'o','MarkerSize',4,'MarkerFaceColor','b')
hold on
plot(tk,dkn,'LineWidth',2)
legend('GPS','Kalman')
grid on
ylabel('Disp. (m)','FontSize',18)
xlim([5 250])
xlabel('Seconds after origin','FontSize',18)
set(gca,'FontSize',16)
subplot(3,3,8)
plot(tgps,egps,'o','MarkerSize',4,'MarkerFaceColor','b')
hold on
plot(tk,dke,'LineWidth',2)
legend('GPS','Kalman')
grid on
xlim([5 250])
xlabel('Seconds after origin','FontSize',18)
set(gca,'FontSize',16)
subplot(3,3,9)
plot(tgps,ugps,'o','MarkerSize',4,'MarkerFaceColor','b')
hold on
plot(tk,dku,'LineWidth',2)
legend('GPS','Kalman')
grid on
xlim([5 250])
set(gca,'FontSize',16)
xlabel('Seconds after origin','FontSize',18)