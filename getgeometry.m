function getgeometry(m,lonmeca,latmeca,zmeca)

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

path='/diego-local/Research/Data/El Mayor'
N=size(m,3);

[m e]=mom2dc(m);
[NP1 NP2 res]=mom2sdr(m);
tplote=1:1:N;
tplot=42:3:300;
tplots=1:1:300;


%Known solutions
%Tohoku
% st1=linspace(193,193,N);
% dip1=linspace(14,14,N);
% rake1=linspace(81,81,N);
% st2=linspace(22,22,N);
% dip2=linspace(76,76,N);
% rake2=linspace(92,92,N);
%Toki
st1=linspace(250,250,300);
dip1=linspace(11,11,300);
rake1=linspace(132,132,300);
st2=linspace(28,28,300);
dip2=linspace(82,82,300);
rake2=linspace(83,83,300);
% %EL Mayor
% st1=linspace(319,319,max(size(tplot)));
% dip1=linspace(77,77,max(size(tplot)));
% rake1=linspace(213,213,max(size(tplot)));
% st2=linspace(222,222,max(size(tplot)));
% dip2=linspace(59,59,max(size(tplot)));
% rake2=linspace(346,346,max(size(tplot)));
% % ______

% %Change sign convention
% i=find(NP2(116:end,1)<180 & NP2(116:end,1)>100)+115;
% NP2(i,1)=NP2(i,1)+180;
% i=find(NP2(116:end,1)<100)+115;
% NP2(i,1)=NP2(i,1)+180;
% i=find(NP1(116:end,1)<150)+115;
% NP1(i,1)=NP1(i,1)+180;
% 
% i=find(NP1(116:end,3)<200)+115;
% NP1(i,3)=NP1(i,3)+180;
% i=find(NP2(116:end,3)<200)+115;
% NP2(i,3)=NP2(i,3)+180;
% 
% %Clean up before threshold
% NP1(1:38,:)=NaN;
% NP2(1:38,:)=NaN;


figure
subplot(4,1,1)
scatter(tplot,NP1(tplot,1),'Marker','v','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
grid on
scatter(tplot,NP2(tplot,1),'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','b')
plot(tplots,st1,'b','LineWidth',2)
plot(tplots,st2,'b','LineWidth',2)
ylabel('Strike (\circ)','FontSize',20)

subplot(4,1,2)
scatter(tplot,NP1(tplot,2),'Marker','v','MarkerFaceColor','gr','MarkerEdgeColor','gr')
hold on
grid on
scatter(tplot,NP2(tplot,2),'Marker','o','MarkerFaceColor','gr','MarkerEdgeColor','gr')
plot(tplots,dip1,'gr','LineWidth',2)
plot(tplots,dip2,'gr','LineWidth',2)
ylabel('Dip \circ','FontSize',20)

subplot(4,1,3)
scatter(tplot,NP1(tplot,3),'Marker','v','MarkerFaceColor','r','MarkerEdgeColor','r')
grid on
hold on
scatter(tplot,NP2(tplot,3),'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r')
plot(tplots,rake1,'r','LineWidth',2)
plot(tplots,rake2,'r','LineWidth',2)
ylabel('Rake \circ','FontSize',20)
%legend('NP1 strike','NP1 dip','NP1 rake','NP2 strike','NP2 dip','NP2 rake')

% figure
% plot(tplot,res)
% grid on
% legend('Nodal plane residual')
subplot(4,1,4)
plot(tplote,e)
grid on
ylabel('\epsilon')
st1=NP1(:,1);
st2=NP2(:,1);
xlabel('Seconds after origin time','FontSize',20)

for k=1:N
     Mo(k)=norm(m(:,:,k),'fro')/sqrt(2); %Scalar moment
     Mw(k)=0.67*(log10(Mo(k))-9.1);
end
writefile=[lonmeca' latmeca' zmeca' Mw' NP1 NP2];
cd(path)
save('fgeometry.dat','writefile','-ascii');
