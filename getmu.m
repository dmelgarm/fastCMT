function mu=getmu(velmod,depth)

%Extract rigidity from velocity model
%
%DMM 05/2011

stop=0;
k=1;
while stop==0
    if depth>velmod(k,1) & depth<velmod(k+1,1)  %interval found
        mu=velmod(k,4)*velmod(k,3)^2;
        stop=1;
    end
    if depth==velmod(k,1)
        mu1=velmod(k,4)*velmod(k,3)^2;
        mu2=velmod(k+1,4)*velmod(k+1,3)^2;
        mu=(mu1+mu2)/2;
        stop=1;
    end
    k=k+1;
end