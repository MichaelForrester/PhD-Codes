% Calculates I'(0) for Jansen-Rit node with parameters P.

function Idash=Idash(P)

% Initial conditions
U0=rand(1,6);
yi=rand(1,6);

% Calculates phase response curve
[Period,~,Y,U] = PRC(P,U0,yi);
u=U(:,5);
clear U;

% Calculates interaction
H=Y(:,2)-Y(:,3);
clear Y;

% Computing derivative
dt=Period/(length(H)-1);
udash=diff(u)/dt;
udash(end+1)=(u(1)-u(end))/dt;
Idash=-udash'*sigm(H)*dt/Period;

end