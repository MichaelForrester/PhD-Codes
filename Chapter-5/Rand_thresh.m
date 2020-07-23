function RandC=Rand_thresh(C,steps)

% Randomise weighted, threholded matrix

for n=1:steps
    
    [I,J]=find(C~=0);
    c1=randi(length(I));
    a=I(c1); b=J(c1);
    cs=find(C(:,b)==0); cs(cs==b)=[];
    ds=find(C(:,a)~=0); ds(ds==b)=[];
    if ~isempty(cs)&&~isempty(ds)
        cdc=C(ds,cs); cdc(cdc==0)=inf;
        cds=find(1-C(ds,a)>cdc|1-C(ds,a)>repmat(C(a,b),size(cdc)));
        if ~isempty(cds)
            c2=cds(randi(length(cds)));
            [di,ci] = ind2sub(size(cdc),c2);
            d=ds(di); c=cs(ci);

            D=min([C(a,b),C(c,d)]);

            C(a,b)=C(a,b)-D;  C(b,a)=C(b,a)-D;
            C(c,d)=C(c,d)-D;  C(d,c)=C(d,c)-D;
            C(a,d)=C(a,d)+D;  C(d,a)=C(d,a)+D;
            C(c,b)=C(c,b)+D;  C(b,c)=C(b,c)+D;
        end
    end
end

RandC=C;