function [t,u] = WC_timeseties(P,C,X0,tspan)

N=length(C); % Number of nodes.

options=[];
[~,X]=ode45(@WCRHS,[0 200],X0,options,P,C); % Initial transients.
X0=X(end,:);
[t,X]=ode45(@WCRHS,tspan,X0,options,P,C); % Time series.
u=X(:,1:N);

end

% Wilson-Cowan ODE system RHS.
function WCRHS=WCRHS(t,X,P,C)

N=length(C); % Number of nodes.

u=X(1:N); v=X(N+1:2*N);

WCRHS=zeros(2*N,1);

WCRHS(1:N)=-u+f(P.c1*u-P.c2*v+P.p+0.01*C*u);
WCRHS(N+1:2*N)=-v+f(P.c3*u-P.c4*v+P.q);

end

% Sigmoid
function f = f(z)

f=1.0./(1+exp(-z));

end