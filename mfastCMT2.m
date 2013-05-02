function [m moment L synthetics lambda_corner]=mfastCMT2(coseis,G,velmod,epi,psmecaf,weightflag,plotflag,dweight,Ln,dcflag,regflag,lambda,tikh)

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
% dcflag - 4 to force Double couple solution
%          5 for full deviatoric MT solution
%wiggle
% Output variables:
%
% m - 3x3xNt array containing the moment tensor at each time step.
% L - 1xNt array containing misfits at each time step.
% synthetics - Observed and synthetic displacements
%_______________________________________________________________________
%_______________________________________________________________________


GMTpath='/diego-local/scripts/GMT/fastCMT'
%Define coordinate origin as current inversion node
late = epi(:,2);
lone = epi(:,1);
depth= epi(:,3);
mu=0;

%Assign data
N=coseis.N';
E=coseis.E';
U=coseis.U';
T=coseis.T';
stdn=coseis.stdn;
stde=coseis.stde;
stdu=coseis.stdu;
lat=coseis.lat;
lon=coseis.lon;

%Get station info
nepi=size(epi,1);
for k=1:nepi
    [dist(:,k),az(:,k)] = distance(late(k),lone(k),lat,lon);
    dist(:,k)=deg2km(dist(:,k))*1000;
    az(:,k)=-az(:,k)+90;
    az(:,k)=deg2rad(az(:,k));
    [e n]=pol2cart(az(:,k),dist(:,k));
    x(:,k)=n;
    y(:,k)=e;
    dsta(:,k)=dist(:,k);
    az(:,k)=atan2(y(:,k),x(:,k));
end
%Get distances between inersion nodes
if nepi>1
    for k=1:nepi-1
        [depi(k),azepi] = distance(late(k),lone(k),late(k+1),lone(k+1));
    end
else
    depi=1;
end
depi=deg2km(depi)*1000;
%Dimensions of data
%nt=322;
nt=size(T,1);
nsta_all=size(T,2);
iall=1:1:nsta_all;
zerod=zeros(size(iall));
ti=1:1:nt;

%Initalize stuff
m=zeros(3,3,nt); %Mwiggleoment tensor
synthetics=zeros(size(dsta,1),8);  %Synthetics stored here
s=1;  %Synthetics matrix counter
thresh=0.0  %Observation threshold




% MAIN PROGRAM _____________________
display('Inverting single CMT...')
mu=getmu(velmod,depth)
dmulti=[];
Gmulti=[];
Wmulti=[];
Wdmulti=[];
for k=1:nt
    for kepi=1:nepi
        Tcurrent=T(k,:);
        Ecurrent=E(k,:);
        Ncurrent=N(k,:);
        Ucurrent=-U(k,:);
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
        
        if nsta>ceil((5/3)*nepi)
            %Prepare data vector
            ds=zeros(nsta*3*nepi,1);
            dx(i)=Ncurrent(i);
            dy(i)=Ecurrent(i);      %Un-normalized for moment
            dz(i)=-Ucurrent(i);
            izero=setxor(i,iall); %Stations to be set to zero because they are under the threshold
            dx(izero)=zerod(izero);
            dy(izero)=zerod(izero);      %Un-normalized for moment
            dz(izero)=zerod(izero);
            i=iall;
            %azi=az(:,kepi);  Don't need this
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
            Gj=G(j,:,kepi);
            %Assemble into column vector
            d(i1)=dz;
            d(i2)=dx;
            d(i3)=dy;
            %Rotate  to z,x,y
%             rot=buildrotmat(azi,2,1);
%             Gj=rot*Gj;
            Gmulti=[Gmulti Gj];
        else %No synthetics
            L(k)=0;
            Mo(k)=0;
            Mw(k)=0;
            synthetics(s,2)=ti(k);
            s=s+1;
        end
    end
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
        Gmulti=W*Gmulti;  %Weigh GFs by same factor as data
    end
    
    if dweight==1   %Weight by distance to source
        kepi=floor(nepi/2);  %Use central node for weighing
        dw=zeros(1,nsta_all*3);
        wd=(dist(i,kepi).^2)/(min(dist(i,kepi))^2);
        dw(i1)=wd; %z
        dw(i2)=wd; %x
        dw(i3)=wd; %y
        Wd=diag(dw);
        d=Wd*d;
        Gmulti=Wd*Gmulti;  %Weigh GFs by same factor as data
    else
        Wd=eye(nsta_all*3);
    end
    %First derivative regularization
    if regflag==1
        R=mderiv(nepi,depi); %First derivative regularization
%         R=lambda*eye(nepi*5);
%         dreg=zeros(nepi*5,1);
%         Gmulti=vertcat(Gmulti,R);
%         d=vertcat(d,dreg);
    end
    
    if Ln==2   %L2
        if dcflag==5
            nsources=size(G,3);
            %Constrain edges to zero
            kz=1e-5;
            Z=zeros(10,nsources*5);
            Z(1,1)=kz;
            Z(2,nsources+1)=kz;
            Z(3,2*nsources+1)=kz;
            Z(4,3*nsources+1)=kz;
            Z(5,4*nsources+1)=kz;
            Z(6,nsources)=kz;
            Z(7,2*nsources)=kz;
            Z(8,3*nsources)=kz;
            Z(9,4*nsources)=kz;
            Z(10,5*nsources)=kz;
            z=zeros(10,1);
            %Set inequalities
            lb=ones(nsources*5,1)*-1e20;
            lb(1,1)=-1e6;
            lb(nsources+1,1)=-1e6;
            lb(2*nsources+1,1)=-1e6;
            lb(3*nsources+1,1)=-1e6;
            lb(4*nsources+1,1)=-1e6;
            lb(nsources,1)=-1e6;
            lb(2*nsources,1)=-1e6;
            lb(3*nsources,1)=-1e6;
            lb(4*nsources,1)=-1e6;
            lb(5*nsources,1)=-1e6;
            ub=-lb;
     
            i1=1:5:nsources*5-4;
            i2=2:5:nsources*5-3;
            i3=3:5:nsources*5-2;
            i4=4:5:nsources*5-1;
            i5=5:5:nsources*5;
            Gi(:,1:nsources)=Gmulti(:,i1);
            Gi(:,nsources+1:2*nsources)=Gmulti(:,i2);
            Gi(:,2*nsources+1:3*nsources)=Gmulti(:,i3);
            Gi(:,3*nsources+1:4*nsources)=Gmulti(:,i4);
            Gi(:,4*nsources+1:5*nsources)=Gmulti(:,i5);

            %Regularization matrix
            %Zeroth order Tikhonov
            if tikh==0
                R=eye(nsources*5);
                dreg=zeros(nsources*5,1);
                [U s_l V]=csvd(Gi);   %Get SVD
            elseif tikh==1%1st order Tikhonov
                R=eye(5*nsources);
                Rd=full(get_l(nsources,1));
                R=blkdiag(Rd,Rd,Rd,Rd,Rd);
                dreg=zeros(size(R,1),1);
                [U s_l V]=cgsvd(Gi,R);   %Get SVD
            elseif tikh==2
                R=eye(5*nsources);
                Rd=full(get_l(nsources,2));
                R=blkdiag(Rd,Rd,Rd,Rd,Rd);
                dreg=zeros(size(R,1),1);
                [U s_l V]=cgsvd(Gi,R);   %Get SVD
            end
            %Iterate
%             kiter=30;
%             lambda=logspace(-9,-13,kiter);
%             %lambda=1e-9;
%             for nreg=1:kiter
%                 %Force inversion
%                 %GR=vertcat(Gi,R*lambda(nreg));
%                 GR=vertcat(Gi,Z,R*lambda(nreg));
%                 %dr=vertcat(d,dreg);
%                 dr=vertcat(d,z,dreg);
%                 [mtinv(:,nreg) resnorm(nreg)]=lsqlin(GR,dr);
%                 %x = lsqlin(GR,dr,[],[],[],[],lb,ub)
%                 dsynth=Gi*mtinv(:,nreg);
%                 %resnorm(nreg)=norm(d-dsynth);
%                 mtnorm(nreg)=norm(R*mtinv(:,nreg));
%             end
%             
%              [lambda_corner,resnorm_corner,mtnorm_corner] = l_corner(resnorm',mtnorm',lambda');
%              lambda_corner=lambda_corner*10
%             loglog(resnorm,mtnorm)
%             hold on
%             scatter(resnorm_corner,mtnorm_corner)
%             GR=vertcat(Gi,Z,R*lambda_corner);
%             [mtinv resnorm]=lsqlin(GR,dr);
%             best=find(resnorm>10);
%             mtnorm=mtnorm(best(end));
%             resnorm=resnorm(best(end))
%             mtinv=mtinv(:,best(end));
% %             figure
            %subplot(2,1,1)
%             loglog(resnorm,mtnorm),grid on
%             ylabel('||Rm||')
%             subplot(2,1,2)
%             loglog(resnorm,lambda)
%             xlabel('||Am-b||')
%             ylabel('\lambda')
%             grid on
%             a=0;


             %[lambda_corner,rho,eta,lambda]=l_curve(U,s_l,d);
             lambda_corner=1e-11;%lambda_corner * 10
             [mtinv,rho,eta]=tikhonov(U,s_l,V,d,lambda_corner);
            
            mt(i1,1)=mtinv(1:nsources,1);
            mt(i2,1)=mtinv(nsources+1:2*nsources,1);
            mt(i3,1)=mtinv(2*nsources+1:3*nsources,1);
            mt(i4,1)=mtinv(3*nsources+1:4*nsources,1);
            mt(i5,1)=mtinv(4*nsources+1:5*nsources,1);
        else
            mt=lsqlin(Gj,d,[],[],[0 0 1 0 0],[0])
        end
    elseif Ln==1   %L1
        if dcflag==5
            m0=ones(5*nepi,1)*1e21/mu;
            mt = l1decode_pd(m0, Gmulti, [], d);
        else
            
        end
    else
        display('FATAL ERROR, no norm selected, must select norm=1 for L1 or norm=2 for L2, terminating execution.')
        return
    end
    % _____________
    
    
    %Compute misfits
    %Rescale for synthetics
    if weightflag==1
        Gmulti=W\Gmulti;
        Gmulti=Wd\Gmulti;
        d=W\d;
        d=Wd\d;
    end
    dsynth=Gmulti*mt;
    L(k)=(1-sqrt((sum((d-dsynth).^2)/(sum(d.^2)))))*100
    %L(k)=norm(dsynth-d)
    
    %Construct synthetics
    %1=station, 2=time, 3,4,5=observed data (Z,X,Y), 6,7,8=inverted
    %data, 9,10, station coords
    Ns=s+nsta_all-1;
    synthetics(s:Ns,1)=i;
    synthetics(s:Ns,2)=ti(k);
    synthetics(s:Ns,3)=d(iinv(1,:));
    synthetics(s:Ns,4)=d(iinv(2,:));
    synthetics(s:Ns,5)=d(iinv(3,:));
    synthetics(s:Ns,6)=dsynth(iinv(1,:));
    synthetics(s:Ns,7)=dsynth(iinv(2,:));
    synthetics(s:Ns,8)=dsynth(iinv(3,:));
    synthetics(s:Ns,9)=lon(i);
    synthetics(s:Ns,10)=lat(i);
    dzobs=d(iinv(1,:));
    dxobs=d(iinv(2,:));
    dyobs=d(iinv(3,:));
    dzsyn=dsynth(iinv(1,:));
    dxsyn=dsynth(iinv(2,:));
    dysyn=dsynth(iinv(3,:));
    %Scale moment tensor by mu and write in cartesian form
    %Project onto line
    lonproj=lone-lone(1);
    latproj=late-late(1);
    s=[lonproj(end);latproj(end)];
    for k=1:nsources
        v=[lonproj(k); latproj(k)];
        proj=(dot(s,v)/dot(s,s))*s;
        moment(k,1)=lone(1)+proj(1);
        moment(k,2)=late(1)+proj(2);
    end
    for kepi=1:nepi
        m(:,:,kepi)=mtinv2mt(mt((kepi-1)*5+1:kepi*5)*mu);
        M(kepi)=norm(m(:,:,kepi),'fro')/sqrt(2); %Scalar moment
        mnorm(:,:,kepi)=m(:,:,kepi)/M(kepi);
        Mw(kepi)=0.67*(log10(M(kepi))-9.1); %Moment magnitude, 9.1 for Nm, 16.1 for dyn-cm
    end
    moment(:,3)=M;%-(0.95*min(M));
    %Make weighted aerage of MTs
    weights=M/sum(M);
    
    ma=mnorm(1,1,:);
    ma=reshape(ma,1,nepi);
    mtavg(1,1)=sum(weights.*ma)/sum(weights);
    
    ma=mnorm(2,2,:);
    ma=reshape(ma,1,nepi);
    mtavg(2,2)=sum(weights.*ma)/sum(weights);
    
    ma=mnorm(3,3,:);
    ma=reshape(ma,1,nepi);
    mtavg(3,3)=sum(weights.*ma)/sum(weights);
    
    ma=mnorm(1,2,:);
    ma=reshape(ma,1,nepi);
    mtavg(1,2)=sum(weights.*ma)/sum(weights);
    mtavg(2,1)=mtavg(1,2);
    
    ma=mnorm(1,3,:);
    ma=reshape(ma,1,nepi);
    mtavg(1,3)=sum(weights.*ma)/sum(weights);
    mtavg(3,1)=mtavg(1,3);
    
    ma=mnorm(2,3,:);
    ma=reshape(ma,1,nepi);
    mtavg(2,3)=sum(weights.*ma)/sum(weights);
    mtavg(3,2)=mtavg(2,3);
    
    mtavg=mtavg/(norm(mtavg,'fro')/sqrt(2));
    m=mtavg*sum(M);
    
%     figure
%     plot(M)% cd(GMTpath);
% save('moment.xyz','moment','-ascii')

%     grid on
%     xlabel('Source #')
%     ylabel('Source Moment (N-m)')
    sum(M)
s=Ns+1;  %Keep track of No of synthetics produced
    
%Save moment function for GMT plotting
% cd(GMTpath);
% save('moment.xyz','moment','-ascii')

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
