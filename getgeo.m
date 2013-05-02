function getgeo(M,lonmeca,latmeca,zmeca)

%DMM 05/2011
% 
% Extract fault geometry information (strike, dip and rake) from momoent tensor
% 
% IN:
% M ~ Array of momeont tensors
% 
% OUT:
% st1 ~ Strike of first nodal plane
% st2 ~ Strike of 2nd nodal plane
%  and a bunch of plots...

close all
path='/diego-local/Research/Data/El Mayor'
N=size(M,3);

[M e]=mom2dc(M);
[NP1 NP2 res]=mom2sdr(M);
%tplot=1:1:N;
tplot=5:5:595;


%Known solutions
%Tohoku
st=linspace(193,193,N);
dip=linspace(14,14,N);
rake=linspace(81,81,N);
%Toki
% st=linspace(250,250,max(size(tplot)));
% dip=linspace(11,11,max(size(tplot)));
% rake=linspace(132,132,max(size(tplot)));
% %EL Mayor
% st=linspace(319,319,max(size(tplot)));
% dip=linspace(77,77,max(size(tplot)));
% rake=linspace(213,213,max(size(tplot)));
% % ______


figure
subplot(3,1,1)
scatter(tplot,NP1(tplot,1),'Marker','v','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
grid on
scatter(tplot,NP2(tplot,1),'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b')
plot(tplot,st,'b','LineWidth',2)
plot(tplot,st-180,'b','LineWidth',2)
ylabel('Strike (\circ)','FontSize',20)

subplot(3,1,2)
scatter(tplot,NP1(tplot,2),'Marker','v','MarkerFaceColor','gr','MarkerEdgeColor','gr')
hold on
grid on
scatter(tplot,NP2(tplot,2),'Marker','o','MarkerFaceColor','gr','MarkerEdgeColor','gr')
plot(tplot,dip,'gr','LineWidth',2)
ylabel('\circ','FontSize',20)

subplot(3,1,3)
scatter(tplot,NP1(tplot,3),'Marker','v','MarkerFaceColor','r','MarkerEdgeColor','r')
grid on
hold on
scatter(tplot,NP2(tplot,3),'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r')
plot(tplot,rake,'r','LineWidth',2)
xlabel('Seconds after origin time','FontSize',20)
ylabel('\circ','FontSize',20)
%legend('NP1 strike','NP1 dip','NP1 rake','NP2 strike','NP2 dip','NP2 rake')

% figure
% plot(tplot,res)
% grid on
% legend('Nodal plane residual')
% figure
% plot(tplot,e)
% grid on
% title('CLVD component')
st1=NP1(:,1);
st2=NP2(:,1);

for k=1:N
     Mo(k)=norm(M(:,:,k),'fro')/sqrt(2); %Scalar moment
     Mw(k)=0.67*(log10(Mo(k))-9.1);
end
writefile=[lonmeca' latmeca' zmeca' Mw' NP1 NP2];
cd(path)
save('fgeometry.dat','writefile','-ascii');
