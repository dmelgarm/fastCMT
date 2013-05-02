function R=mderiv(nsources,dist)

nparams=5;
R=[];
dist=dist;
dist(:)=1;
for k=0:nsources-1;
    Rd=zeros(nparams,nparams*nsources);
    if k<nsources-1
        for j=1:nparams
            i1=j+k*nparams;
            i2=j+(k+1)*nparams;
            Rd(j,i1)=-1/dist(k+1);
            Rd(j,i2)=1/dist(k+1);
        end
    else
        for j=1:nparams
            i1=j+(k-1)*nparams;
            i2=j+k*nparams;
            Rd(j,i1)=1/dist(k);
            Rd(j,i2)=-1/dist(k);
        end
    end
    R=vertcat(R,Rd);
end