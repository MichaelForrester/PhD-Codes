function RandC=RandC(C,steps)

% Randomise weighted, all-to-all matrix.

N=length(C);

I=zeros(steps,4);
I(:,1)=randi(N,steps,1);

I(:,2)=randi(N-1,steps,1); I(I(:,2)==I(:,1),2)=N;

I(:,3)=randi(N-2,steps,1);
I(I(:,3)>=min(I(:,1:2),[],2),3)=I(I(:,3)>=min(I(:,1:2),[],2),3)+1;
I(I(:,3)>=max(I(:,1:2),[],2),3)=I(I(:,3)>=max(I(:,1:2),[],2),3)+1;

I(:,4)=randi(N-3,steps,1);
I(I(:,4)>=min(I(:,1:3),[],2),4)=I(I(:,4)>=min(I(:,1:3),[],2),4)+1;
I(I(:,4)>=median(I(:,1:3),2),4)=I(I(:,4)>=median(I(:,1:3),2),4)+1;
I(I(:,4)>=max(I(:,1:3),[],2),4)=I(I(:,4)>=max(I(:,1:3),[],2),4)+1;

for n=1:steps
    a=I(n,1); b=I(n,2); c=I(n,3); d=I(n,4);
    
    D=min([C(a,b),C(c,d),1-C(a,d),1-C(b,c)]);
    
    C(a,b)=C(a,b)-D;  C(b,a)=C(b,a)-D;
    C(c,d)=C(c,d)-D;  C(d,c)=C(d,c)-D;
    C(a,d)=C(a,d)+D;  C(d,a)=C(d,a)+D;
    C(c,b)=C(c,b)+D;  C(b,c)=C(b,c)+D;   
end

RandC=C;