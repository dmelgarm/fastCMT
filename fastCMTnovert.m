function [m L synthetics]=fastCMTnovert(coseis,G,velmod,epi,psmecaf,weightflag,plotflag,dweight,Ln)

%DMM 01/2011
%
% Invert for moment tensor each time sample using GPS data. Inversion is
% performed using the coseismic offsets with a Green's function approach.
% Green's functions are obtained from EDGRN numerical code (Wang et al.,
% 2003).
%
% Coordinate system throughout this code and the functions called within is
% (X,Y,Z)=(North,East,Down), this is important, save yourself headaches and 
% remember this fact.
%
% Input variables:
%
% coseis - Structure variable containing the following fields
%     .N ~ North dispalcement;
%     .E ~ East dispalcement;
%     .U ~ Vertical dispalcement;
%     .T ~ Time, zero is assumed to be the event origin time;
%     .lon ~ Station longitudes;
%     .lat ~ Station latitudes;
%     .stdn ~ Uncertainties in north observations (pre-event std. deviaiton);
%     .stde ~ Uncertainties in east observations;
%     .stdu ~ Uncertainties in vertical observations;
% G - Green function matrix for all nodes
% velmod - Four column strucure model, column 1 is depth of layer, column 2
%          is p-wave velocity, column 3 is s wave velocity, column 4 is density
%          like in many other codes, boundaries between layer are repeated.
% epi - epicenter [lon lat depth(m)]
% psmecaf - psmeca filename
% weightflag - Weight by preevent standard deviations on all 3 channels.
% plotflag - Plot observed andsynthetic offsets each time sample.
% dweight - Weight by 1/d^2 (distance froms tation to inversion node
% Ln - 1 for L1 inversion, 2 for L2 inversion
%
% Output variables:
%
% m - 3x3xNt array containing the moment tensor at each time step.
% L - 1xNt array containing misfits at each time step.
% synthetics - Observed and synthetic displacements
%_______________________________________________________________________
%_______________________________________________________________________



%Define coordinate origin as current inversion node
late = epi(2);
lone = epi(1);
depth= epi(3);
mu=0;

%Assign data
N=coseis.N';
E=coseis.E';
U=-coseis.U';
T=coseis.T';
stdn=coseis.stdn;
stde=coseis.stde;
stdu=coseis.stdu;
lat=coseis.lat;
lon=coseis.lon;

%Get station info
[dist,az] = distance(late,lone,lat,lon);
dist=deg2km(dist)*1000;
az=-az+90;
az=deg2rad(az);
[e n]=pol2cart(az,dist);
x=n;
y=e;
dsta=dist;
az=atan2(y,x);

%Dimensions of data
%nt=322;
nt=size(T,1);
nsta_all=size(T,2);
iall=1:1:nsta_all;
zerod=zeros(size(iall));
ti=1:1:nt;

%Initalize stuff
m=zeros(3,3,nt); %Moment tensor
synthetics=zeros(size(dsta,1),8);  %Synthetics stored here
s=1;  %Synthetics matrix counter
thresh=0.0  %Observation threshold




% MAIN PROGRAM _____________________
display('Inverting single CMT...')
mu=getmu(velmod,depth)
for k=1:nt
    Tcurrent=T(k,:);
    Ecurrent=E(k,:);
    Ncurrent=N(k,:);
    Ucurrent=U(k,:);
    i=find(~isnan(Ncurrent));
    %Use only stations over threshold
    horiz=sqrt(Ecurrent.^2+Ncurrent.^2);
    i=find(horiz>=thresh);
    nsta=max(size(i));
    %Fix stations after max disp reached
%     if k<600
%         i=find(horiz>=thresh);
%         if ~isempty(i)
%             nsta=max(size(i));  %Stations in this epoch
%         else
%             nsta=0;
%         end
%     elseif k==600
%         imax=find(horiz>=thresh);
%         i=imax;
%         if ~isempty(i)
%             nsta=max(size(i));  %Stations in this epoch
%         else
%             nsta=0;
%         end
%     else
%         i=imax;
%     end
    %Noise matrix
    W=eye(nsta*3);
    if nsta>2
        %Prepare data vector
        ds=zeros(nsta*3,1);
        dx(i)=Ncurrent(i);
        dy(i)=Ecurrent(i);      %Un-normalized for moment
        dz(i)=Ucurrent(i);
        izero=setxor(i,iall); %Stations to be set to zero because they are under the threshold
        dx(izero)=zerod(izero);
        dy(izero)=zerod(izero);      %Un-normalized for moment
        dz(izero)=zerod(izero);
        i=iall;
        azi=az(i);
        d=ones(size(i,1)*3,1);
        iinv=1:1:nsta_all*3;
        iinv=reshape(iinv,3,nsta_all);
        %Indices for preparing d vector
        i1=iinv(1,:);
        i2=iinv(2,:);
        i3=iinv(3,:);
        %Extract only relevant rows of Gram matrix
        ii1=i*3-2;
        ii2=i*3-1;
        ii3=i*3;
        j1=interleave(ii1,ii2);
        j2=interleave(ii2,ii3);
        %j=unique(interleave(j1,j2));
        j=j2;
        Gj=G(j,:);
        %redefine counters for no vert
        iinv=1:1:nsta_all*2;
        iinv=reshape(iinv,2,nsta_all);
        i1=iinv(1,:);
        i2=iinv(2,:);
        %Assemble data and get rid of verticals
        d(i1)=dx;
        d(i2)=dy;
        %Rotate un-normalized data to z,r,theta
        rot=buildrotmat(azi,2,0);
        d=rot*d;
        %Weight by standard errors
        if weightflag==1
            %Simplified weight scheme
            dw=zeros(1,nsta_all*2);
            %dw(i1)=1./stdu(i); %z
            dw(i1)=1./stdn(i); %x
            dw(i2)=1./stde(i); %y
            W=diag(dw);
            is=isnan(W);
            W(is)=0;
            d=W*d;
            Gj=W*Gj;  %Weigh GFs by same factor as data
        end
        
        if dweight==1   %Weight by distance to source
            dw=zeros(1,nsta_all*2);
            wd=(dist(i).^2)/(min(dist(i))^2);
            dw(i1)=wd; %x
            dw(i2)=wd; %y
            Wd=diag(dw);
            d=Wd*d;
            Gj=Wd*Gj;  %Weigh GFs by same factor as data
        else
            Wd=eye(nsta_all*2);
        end
            
        if Ln==2   %L2
            mt=lsqlin(Gj,d);
        elseif Ln==1   %L1
            m0=ones(5,1)*1e17/mu;
            mt = l1decode_pd(m0, Gj, [], d);
        else
            display('FATAL ERROR, no norm selected, must select norm=1 for L1 or norm=2 for L2, terminating execution.')
            return
        end
        % _____________

        

        %Compute misfits
        %Rescale for synthetics
        if weightflag==1
            Gj=W\Gj;
            Gj=Wd\Gj;
            d=W\d;
            d=Wd\d;
        end
        dobs=d;
        dsyn=Gj*mt; 
        L(k)=(1-(sqrt(sum((dobs-dsyn).^2)/(sum(dobs.^2)))))*100;
        
        %Synthetics
        dsynth=Gj*mt;
        %Rotate back to z,x,y
        rot=buildrotmat(azi,1,0);
        dsynth=rot*dsynth;
        d=rot*d;

        %Construct synthetics
        %1=station, 2=time, 3,4,5=observed data (Z,X,Y), 6,7,8=inverted
        %data, 9,10, station coords
        Ns=s+nsta_all-1;
        synthetics(s:Ns,1)=i;
        synthetics(s:Ns,2)=ti(k);
        synthetics(s:Ns,3)=zeros(size(d(iinv(1,:))));
        synthetics(s:Ns,4)=d(iinv(1,:));
        synthetics(s:Ns,5)=d(iinv(2,:));
        synthetics(s:Ns,6)=zeros(size(dsynth(iinv(1,:))));
        synthetics(s:Ns,7)=dsynth(iinv(1,:));
        synthetics(s:Ns,8)=dsynth(iinv(2,:));
        synthetics(s:Ns,9)=lon(i);
        synthetics(s:Ns,10)=lat(i);
        %Scale moment tensor by mu and write in cartesian form
        m(:,:,k)=mtinv2mt(mt*mu);
        M(k)=norm(m(:,:,k),'fro')/sqrt(2); %Scalar moment
        Mw(k)=0.67*(log10(M(k))-9.1); %Moment magnitude, 9.1 for Nm, 16.1 for dyn-cm
        s=Ns+1;  %Keep track of No of synthetics produced
    else %No synthetics
        L(k)=0;
        Mo(k)=0;
        Mw(k)=0;
        synthetics(s,2)=ti(k);
        s=s+1;
    end
end

if plotflag==1
    clear G N E U T
    %Make plots of syntehtics
    figure
    tplot=unique(synthetics(:,2));
    for k=1:max(size(tplot,1))
        
        [row col]=find(synthetics(:,2)==tplot(k));
        t(k)=tplot(k);
        if synthetics(row(1),1)>0
            
            ista=synthetics(row,1);
            dz=synthetics(row,3);
            dx=synthetics(row,4);
            dy=synthetics(row,5);
            sz=synthetics(row,6);
            sx=synthetics(row,7);
            sy=synthetics(row,8);
            
            
            %Get station coordinates
            xc=x(ista)/1000;
            yc=y(ista)/1000;
            
            %Plot, plot, plotaroo...
            scale=2;
            subplot(1,2,1)   %Horizontals
            quiver(yc,xc,dy,dx,scale)
            hold on
            quiver(yc,xc,sy,sx,scale)
            plot(0,0,'x','MarkerSize',20)
            
            %Plot limits
            limx=[min(yc) max(yc)];
            dlimx=(limx(2)-limx(1))*0.1;
            limx(1)=limx(1)-dlimx;
            limx(2)=limx(2)+dlimx;
            limy=[min(xc) max(xc)];
            dlimy=(limy(2)-limy(1))*0.1;
            limy(1)=limy(1)-dlimy;
            limy(2)=limy(2)+dlimy;
            
            xlim(limx);
            ylim(limy);
            axis equal
            grid on
            xlabel('East (km)','FontSize',18)
            ylabel('North (km)','FontSize',18)
            title(['t = ' num2str(t(k))],'FontSize',18);
            legend('Observed','Synthetic')
            hold off
            subplot(1,2,2)   %verticals
            pq=zeros(size(xc));
            quiver(yc,xc,pq,dz,scale)
            hold on
            quiver(yc+5,xc,pq,sz,scale)
            plot(0,0,'x','MarkerSize',20)
            
            xlim(limx);
            ylim(limy);
            axis equal
            grid on
            xlabel('East (km)','FontSize',18)
            ylabel('North (km)','FontSize',18)
            title(['t = ' num2str(t(k))],'FontSize',18);
            legend('Observed','Synthetic')
            hold off
            
            pause(0.01)
            
            
        else
            subplot(1,2,1)
            xlim([min(x) max(x)]);
            ylim([min(y) max(y)]);
            axis equal
            grid on
            xlabel('East (km)','FontSize',18)
            ylabel('North (km)','FontSize',18)
            title(['t = ' num2str(t(k))],'FontSize',18);
            subplot(1,2,2)
            xlim([min(x)/1000 max(x)/1000]);
            ylim([min(y)/1000 max(y)/1000]);
            axis equal
            grid on
            xlabel('East (km)','FontSize',18)
            ylabel('North (km)','FontSize',18)
            title(['t = ' num2str(t(k))],'FontSize',18);
            pause(0.01)
        end
    end
end
