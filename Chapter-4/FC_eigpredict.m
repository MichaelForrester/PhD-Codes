parpool('local',40)

tic;
load('Connect2.mat')
Connect_bin=Connect;
Connect=Connect_bin./sum(Connect_bin,2);
F=size(Connect,1);

maxruns=1000;

P.C1=135;
P.C2=0.8*135;
P.C3=0.25*135;
P.C4=0.25*135;
P.a=100; P.b=50;
P.P=120;

A=6; B=18;

ep=1;

stab=zeros(101);

T = 100; dt=0.0001; tspan=0:dt:T;

U0=rand(1,6);
yi=rand(1,6);
[Period,~,Y,U] = PRC_JRnew(A,B,U0,yi);
u=U(:,5);
U=[];
IC=zeros(maxruns,6*F);
for n=1:maxruns
    y_i=Y(round(size(Y,1)*rand(1,F)),:);
    IC(n,:)=y_i(:);
end
y_prc=Y(:,2)-Y(:,3);
Y=[];

FC_sim=zeros(F,F,maxruns);
FC_pred=zeros(F,F,maxruns);

parfor runs=1:maxruns
    ie=2;
    while length(ie)<F
        [~,y]=ode45(@(t,y) JansenRit(t,y,A,B,P,ep,Connect,F), 0:0.001:T, IC(n,:));
        y_i=y(end,:);
        opts = odeset('Events',@(t,X) PDFind(t,X,y_i(1:F)),'RelTol',1e-6,'AbsTol',1e-8);
        [~,~,te,~,ie] =ode45(@(t,y) JansenRit(t,y,A,B,P,ep,Connect,F), [0 5*Period], y_i,opts);
    end
    y0save=y(:,1:F);
    y=[];

    U_trans = imag(hilbert(y0save));
    R = zeros(F);
    N=length(tspan);
    for f = 1:F
                R(f,f+1:F) = abs((1/N)*sum(exp(1i*(U_trans(:,f+1:F)-repmat(U_trans(:,f),1,F-f)))));
                R(f+1:F,f)=R(f,f+1:F);
    end
    
    FC_sim(:,:,runs)=R;
    
    PD=PDForm(ie,te,F,Period);
    [V,D]=eig(H_Jac_New(PD,A,B,y_prc,u,Period));
    [~,I] = min(abs(real(diag(D))));
    V(:,I)=[]; D(:,I)=[];
    [~,I] = max(abs(real(diag(D))));
    FC_pred(:,:,runs)=V(:,I)*V(:,I)';
end

FC_pred=sum(FC_pred,3)/maxruns;
FC_sim=sum(FC_sim,3)/maxruns;

save('FCs','FC_pred','FC_sim');

toc

poolobj=gcp('nocreate');
delete(poolobj);