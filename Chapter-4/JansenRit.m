% Jansen-Rit ODE

function dydt=JansenRit(t,y,P,C,N)
    
    y0=y(1:N);
    y1=y(N+1:2*N);
    y2=y(2*N+1:3*N);
    y3=y(3*N+1:4*N);
    y4=y(4*N+1:5*N);
    y5=y(5*N+1:6*N);
    
    dydt=zeros(6*N,1);
    
    dydt(1:N)=y3;
    dydt(N+1:2*N)=y4;
    dydt(2*N+1:3*N)=y5;
    dydt(3*N+1:4*N)=P.A*P.a*sigm(y1-y2)-2*P.a*y3-(P.a^2)*y0;
    dydt(4*N+1:5*N)=P.A*P.a*(P.P+P.ep*C*sigm((y1-y2))+P.C2*sigm(P.C1*y0))-2*P.a*y4-(P.a^2)*y1;
    dydt(5*N+1:6*N)=P.B*P.b*P.C4*sigm(P.C3*y0)-2*P.b*y5-(P.b^2)*y2;
    
end