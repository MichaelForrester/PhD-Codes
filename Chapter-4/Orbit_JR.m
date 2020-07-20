%Orbit_WC
function [Period,t,Y] = Orbit(P,yi)

% Total time
total=200;

% Initial data
X0=yi;

% Evolve for large time
[t,X]=ode45(@JRRHS,[0 total],X0,options,P,A,B);

% Use findpeaks to estimate period
[~,b] = findpeaks(X(:,1));

% Calculate period (if it exists)
if ~isempty(b) && length(b) >=20     % Find 20 crossings and then average
    Period = (t(b(end))-t(b(end-20)))/20;
else
    Period = 0;
    t=[];
    Y=[];
    return
end

% Terminates if 

xmax=max(X(ceil(end/2):end,1));
xmin=min(X(ceil(end/2):end,1));

if Period < 0.01 || xmax-xmin < 0.01
    Period = 0;
    t=[];
    Y=[];
    return
end

%run over one mutiple of first guess for period to get better estimate
total=30*Period;
Y0=X(end,:);
%options=['RelTol',1e-6,'AbsTol',1e-7];
options=odeset('Events',@(t,X) PeriodFind(t,X,Y0),'RelTol',1e-8,'AbsTol',1e-10);
[~,X,te,~,~]=ode45(@(t,X) JRRHS(t,X,P,A,B),[0 total],Y0,options);
while length(te)<2
    disp('orbit')
    options=odeset('Events',@(t,X) PeriodFind2(t,X,Y0),'RelTol',1e-8,'AbsTol',1e-10);
    Y0=X(round((end/2)*(1+rand)),:);
    [~,X,te,~,~]=ode45(@(t,X) JRRHS(t,X,P,A,B),[0 total],Y0,options);
end
Period=te(end)-te(end-1);

%one more time for luck! ... really to get even mesh
N=1000000;
dt=Period/N;
total=Period;
Y0=X(end,:);
options=odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,X]=ode45(@(t,X) JRRHS(t,X,P,A,B),0:dt:total,Y0,options);

Y=X;


% 
% figure(10)
% clf
% plot(t,E,'ro-')
% axis([0 Period 0 1])
% 
% figure(11)
% clf
% plot(E,I,'bo-')

end

function JRRHS=JRRHS(t,X,P,A,B)

y=zeros(6,1);

y(1)=X(1);
y(2)=X(2);
y(3)=X(3);
y(4)=X(4);
y(5)=X(5);
y(6)=X(6);

JRRHS=[y(4);y(5);y(6);
    A*P.a*sigm(y(2)-y(3))-2*P.a*y(4)-(P.a^2)*y(1);
    A*P.a*(P.P+P.C2*sigm(P.C1*y(1)))-2*P.a*y(5)-(P.a^2)*y(2);
    B*P.b*P.C4*sigm(P.C3*y(1))-2*P.b*y(6)-(P.b^2)*y(3)];
end

function sigm = sigm(v)

v0=6.0;
e0=2.5;
r=0.56;

sigm=(2*e0)./(1+exp(r*(v0-v)));

end

function [position,isterminal,direction] = PeriodFind(t,X,Y0)
position = X(1)-Y0(1); % The value that we want to be zero
isterminal = 0;  % Halt integration 
direction = 1;   % The zero can be approached from either direction
end

function [position,isterminal,direction] = PeriodFind2(t,X,Y0)
position = X(1)-Y0(1); % The value that we want to be zero
isterminal = 0;  % Halt integration 
direction = -1;   % The zero can be approached from either direction
end