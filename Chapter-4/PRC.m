function [Period,t,Y,U] = PRC(P,U0,yi)

% Construct orbit
[Period,t,Y] = Orbit(P,yi);

% Check for existence of orbit
if Period == 0
    U=[];
    return
else

% Evolve adjoint equation in backward time
total=20*Period;
dt=t(2)-t(1);
options=[];
X0=U0;
[t,X]=ode45(@Adjoint,total:-dt:0,X0,options,P,Y,A,B,Period);

% Normalise over one period
t=t(end-length(Y(:,1))+1:end);
X=X(end-length(Y(:,1))+1:end,:);

U=X;

U=flip(U); % Backwards time

%shift origin of time to zero
t=t-(total-Period);

%compute inner product of adjoint and vector field
ydot=[Y(:,4),Y(:,5),Y(:,6),A*P.a*sigm(Y(:,2)-Y(:,3))-2*P.a*Y(:,4)-(P.a^2)*Y(:,1),A*P.a*(P.P+P.C2*sigm(P.C1*Y(:,1)))-2*P.a*Y(:,5)-(P.a^2)*Y(:,2),B*P.b*P.C4*sigm(P.C3*Y(:,1))-2*P.b*Y(:,6)-(P.b^2)*Y(:,3)];
inner=sum(diag(U'*ydot))*dt/Period;
U=U/inner;

end

end


function Adjoint=Adjoint(t,X,P,Y,A,B,Period)

%backward time: t -> -t
U.u0=X(1);
U.u1=X(2);
U.u2=X(3);
U.u3=X(4);
U.u4=X(5);
U.u5=X(6);

% Number of time steps in the orbit
N=length(Y(:,1));

% Backward time
t=mod(t,Period);
index=floor(t/Period*N)+1;
y0back=Y(index,1);
y1back=Y(index,2);
y2back=Y(index,3);

% Adjoint ODE system
Adjoint=-[-(P.a^2)*U.u3+A*P.a*P.C2*P.C1*dsigm(P.C1*y0back)*U.u4+B*P.b*P.C4*P.C3*dsigm(P.C3*y0back)*U.u5;
    A*P.a*dsigm(y1back-y2back)*U.u3-(P.a^2)*U.u4;
    -A*P.a*dsigm(y1back-y2back)*U.u3-(P.b^2)*U.u5;
    U.u0-2*P.a*U.u3; U.u1-2*P.a*U.u4; U.u2-2*P.b*U.u5];

end

function sigm = sigm(v)

% Sigmoid

v0=6.0;
e0=2.5;
r=0.56;

sigm=(2*e0)./(1+exp(r*(v0-v)));

end

function dsigm = dsigm(v)

% 1st derivative of sigmoid.

e0=2.5;
r=0.56;

dsigm = r*(1-(1/(2*e0))*sigm(v)).*sigm(v);

end