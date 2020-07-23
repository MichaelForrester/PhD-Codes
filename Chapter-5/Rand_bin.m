function RandC=Rand_bin(C,steps)

% Randomise binarised, thresholded matrix

for n=1:steps
    
    [I,J]=find(C~=0);
    c1=randi(length(I));
    a=I(c1); b=J(c1);
    cs=find(C(:,b)==0); cs(cs==b)=[];
    ds=find(C(:,a)==0); ds(ds==a)=[];
    if ~isempty(cs)&&~isempty(ds)
        cdc=C(ds,cs);
        cds=find(cdc==1);
        if ~isempty(cds)
            c2=cds(randi(length(cds)));
            [di,ci] = ind2sub(size(cdc),c2);
            d=ds(di); c=cs(ci);

            C(a,b)=0;  C(b,a)=0;
            C(c,d)=0;  C(d,c)=0;
            C(a,d)=1;  C(d,a)=1;
            C(c,b)=1;  C(b,c)=1;
        end
    end
end

RandC=C;