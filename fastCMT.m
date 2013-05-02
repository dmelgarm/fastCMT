function [m moment L synthetics]=fastCMT(coseis,G,velmod,epi,psmecaf,weightflag,plotflag,dweight,Ln,dcflag,secondflag)

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
%          is p-wave velocity, column 3 is s wave velocity, column 4 iswiggle density
%          like in many other codes, boundaries between layer are repeated.
% epi - epicenter [lon lat depth(m)]
% psmecaf - psmeca filename
% weightflag - Weight by preevent standard deviations on all 3 channels.
% plotflag - Plot observed andsynthetic offsets each time sample.
% dweight - Weight by 1/d^2 (distance froms tation to inversion node
% Ln - 1 for L1 inversion, 2 for L2 inversion
% dcflag - 4 to force Double couple solution
%          5 for full deviatoric MT solution
% secondfalg - 1 for 2nd order MT ivnersion
%              0 for zeroth order (traditional) MT inversion
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
    %Use only stations over threshold
    horiz=sqrt(Ecurrent.^2+Ncurrent.^2);
    i=find(horiz>=thresh);
%     Ecurrent(i)=0;
%     Ncurrent(i)=0;
%     Ucurrent(i)=0;
%     i=find(horiz>=thresh);
    nsta=max(size(i));

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
        j=unique(interleave(j1,j2));
        Gj=G(j,:,k,1);  %Load GFs
        %Assemble into column vector
        d(i1)=dz;
        d(i2)=dx;
        d(i3)=dy;
        du=d; %Unweighted data
        %Rotate GFs to x,y,z
%         rot=buildrotmat(azi,1,1);
%         Gj=rot*Gj;
        Gu=Gj;  %Unweighted GFs
        %Weight by standard errors
        if weightflag==1
            %Simplified weight scheme
            dw=zeros(1,nsta_all*3);
            stdu=stdu*1;
            dw(i1)=1./stdu(i); %z
            dw(i2)=1./stdn(i); %x
            dw(i3)=1./stde(i); %y
            W=diag(dw);
            is=isnan(W);
            W(is)=0;
            d=W*d;
            Gj=W*Gj;  %Weigh GFs by same factor as data
        else
            W=eye(nsta_all*3);
        end
        
        if dweight==1   %Weight by distance to source
            dw=zeros(1,nsta_all*3);
            wd=(dist(i).^2)/(min(dist(i))^2);
            dw(i1)=wd; %z
            dw(i2)=wd; %x
            dw(i3)=wd; %y
            Wd=diag(dw);
            d=Wd*d;
            Gj=Wd*Gj;  %Weigh GFs by same factor as data
        else
            Wd=eye(nsta_all*3);
        end
           
        % Zeroth order MT inversion
        if Ln==2   %L2
            if dcflag==5
                mt=lsqlin(Gj,d);
            else
                mt=lsqlin(Gj,d,[],[],[0 0 1 0 0],[0])
            end
        elseif Ln==1   %L1
            if dcflag==5
                m0=ones(5,1)*1e22/mu;
                mt = l1decode_pd(m0, Gj, [], d);
            else
             
            end
        else
            display('FATAL ERROR, no norm selected, must select norm=1 for L1 or norm=2 for L2, terminating execution.')
            return
        end
        %Compute misfits
        %Rescale for synthetics
        if weightflag==1
            Gj=W\Gj;
            Gj=Wd\Gj;
            d=W\d;
            d=Wd\d;
        end
        d0=Gu*mt; 
        L(k)=(1-((sum((du-d0).^2)/(sum(du.^2)))))*100;
        if secondflag==1
            %Form 2nd order MT GFs
            %Centroid GFs
            G2(:,1)=G(:,:,k,2)*mt;  % x-coordinate of centroid
            G2(:,2)=G(:,:,k,3)*mt;  % y-coordinate of centroid
            G2(:,3)=G(:,:,k,4)*mt;  % z-coordinate of centroid
            %Spatial extent (mu)
            G2(:,4)=G(:,:,k,5)*mt;  % mu_xx
            G2(:,5)=G(:,:,k,6)*mt;  % mu_yy
            G2(:,6)=G(:,:,k,7)*mt;  % mu_zz
            G2(:,7)=G(:,:,k,8)*mt;  % mu_xy
            G2(:,8)=G(:,:,k,9)*mt;  % mu_xz
            G2(:,9)=G(:,:,k,10)*mt; % mu_yz
            
            Gx=G(:,:,k,2);  % x-coordinate of centroid
            Gy=G(:,:,k,3);  % y-coordinate of centroid
            Gz=G(:,:,k,4);  % z-coordinate of centroid
            %Spatial extent (mu)
            Gxx=G(:,:,k,5);  % mu_xx
            Gyy=G(:,:,k,6);  % mu_yy
            Gzz=G(:,:,k,7);  % mu_zz
            Gxy=G(:,:,k,8);  % mu_xy
            Gxz=G(:,:,k,9);  % mu_xz
            Gyz=G(:,:,k,10); % mu_yz
            x=sdpvar(1);
            y=sdpvar(1);
            z=sdpvar(1);
            m=sdpvar(5,1);
            muxx=sdpvar(1);
            muyy=sdpvar(1);
            muzz=sdpvar(1);
            muxy=sdpvar(1);
            muxz=sdpvar(1);
            muyz=sdpvar(1);
            M=[muxx muxy muxz;muxy muyy muyz;muxz muyz muzz];
            objective=norm(du-(Gu+Gx*x+Gy*y+Gz*z+0.5*Gxx*muxx+0.5*Gyy*muyy+0.5*Gzz*muzz+Gxy*muxy+Gxz*muxz+Gyz*muyz)*mt);
            constraints=[M > 0];
            
            e=inf;  %looping tolerance
            tol=1e-4;
            k2=0;
            while e>tol
                k2=k2+1;

                    
                %Second order MT inversion
                %Invert
                display('L2 inversion of 2nd order MT...') 
                %Set up SDP solver
                M2 = sdpvar(3);
                centroid = sdpvar(3,1);
                mt2 = [centroid;M2(find(triu(ones(3))))]
                Objective = norm(G2*mt2-(du-d0));
                Constraints = [M2 > 0, mt2(3)>-depth];
                solvesdp(Constraints,Objective)
                %Convert back to double
                mt2=double(mt2);
                d2=G2*mt2;
                dcorr=du-d2;
                %clean up
                clear centroid M2 Objective Constrints
                %Invert the data corrected for 2nd order MT effects
%                 if weightflag==1  %distance weight
%                     dcorr=W*dcorr;
%                     Gj=W*Gj;  %Weigh GFs by same factor as data
%                 end
%                 if dweight==1   %Weight by distance to source
%                     dcorr=Wd*dcorr;
%                     Gj=Wd*Gj;  %Weigh GFs by same factor as data
%                 else
%                     Wd=eye(nsta_all*3);
%                 end
                dxobs=du(iinv(2,:));
                dyobs=du(iinv(3,:));
                dx2=dcorr(iinv(2,:));
                dy2=dcorr(iinv(3,:));
                dxcorr=d2(iinv(2,:));
                dycorr=d2(iinv(3,:));
                mtemp=mtinv2mt(mt*mu);
                Mtemp=norm(mtemp,'fro')/sqrt(2)
                Mwtemp=0.67*(log10(Mtemp)-9.1)
                mt=lsqlin(W*Gu,W*dcorr); %invert corrected and weighed data
                mtemp=mtinv2mt(mt*mu);
                Mtemp=norm(mtemp,'fro')/sqrt(2)
                Mwtemp=0.67*(log10(Mtemp)-9.1)
                %Re-scale
%                 if weightflag==1
%                     Gj=W\Gj;
%                 end
%                 if dweight==1
%                     Gj=Wd\Gj;
%                 end
                d0corr=Gu*mt;  %Corrected synthetics
                e(k2)=sqrt((d0corr-d0)'*(d0corr-d0));
                d0=d0corr;
                L(k)=(1-((sum((du-d0).^2)/(sum(du.^2)))))*100;
                m(:,:,k)=mtinv2mt(mt*mu);
                M(k)=norm(m(:,:,k),'fro')/sqrt(2); %Scalar moment
                Mw(k)=0.67*(log10(M(k))-9.1); %Moment magnitude, 9.1 for Nm, 16.1 for dyn-cm
%                 %Debugging stuff
%                 dxsyn=d0(iinv(2,:));
%                 dysyn=d0(iinv(3,:));
%                 dxobs=d(iinv(2,:));
%                 dyobs=d(iinv(3,:));
%                 figure
%                 quiver(y,x,dyobs,dxobs,3),hold on,grid on,quiver(y,x,dysyn,dxsyn,3),legend('Obs','Syn'),axis equal   %For debugging
                if k2==10
                    e=1e-10
                    mt2
                end
            end
        end
        %Construct synthetics
        %1=station, 2=time, 3,4,5=observed data (Z,X,Y), 6,7,8=inverted
        %data, 9,10, station coords
        Ns=s+nsta_all-1;
        synthetics(s:Ns,1)=i;
        synthetics(s:Ns,2)=ti(k);
        synthetics(s:Ns,3)=du(iinv(1,:));
        synthetics(s:Ns,4)=du(iinv(2,:));
        synthetics(s:Ns,5)=du(iinv(3,:));
        synthetics(s:Ns,6)=d0(iinv(1,:));
        synthetics(s:Ns,7)=d0(iinv(2,:));
        synthetics(s:Ns,8)=d0(iinv(3,:));
        synthetics(s:Ns,9)=lon(i);
        synthetics(s:Ns,10)=lat(i);
        %Scale moment tensor by mu and write in cartesian form
        m(:,:,k)=mtinv2mt(mt*mu);
        M(k)=norm(m(:,:,k),'fro')/sqrt(2); %Scalar moment
        moment=M;
        Mw(k)=0.67*(log10(M(k))-9.1); %Moment magnitude, 9.1 for Nm, 16.1 for dyn-cm
        s=Ns+1;  %Keep track of No of synthetics produced
        
        %Plot for debugging
%         figure
%         quiver(synthetics(:,9),synthetics(:,10),synthetics(:,5),synthetics(:,4),3)
%         hold on,grid on
%         quiver(synthetics(:,9),synthetics(:,10),synthetics(:,8),synthetics(:,7),3)
%         legend('Obs','Syn')
%         scatter(lone,late,'x')
%         axis equal
        
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
