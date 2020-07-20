function [Period,t,Y] = Orbit(P,yi)

% Total time.
total=200;

% Initial data.
X0=yi;

% Evolve for large time.
[t,X]=ode45(@JR,[0 total],X0,options,P,A,B);

% Use findpeaks to estimate period.
[~,b] = findpeaks(X(:,1));

% Calculate period (if it exists).
if ~isempty(b) && length(b) >=20     % Find 20 crossings and then average.
    Period = (t(b(end))-t(b(end-20)))/20;
else
    Period = 0;
    t=[];
    Y=[];
    return
end

% Terminates if period is too small or solution does not converge.
xmax=max(X(ceil(end/2):end,1));
xmin=min(X(ceil(end/2):end,1));
if Period < 0.01 || xmax-xmin < 0.01
    Period = 0;
    t=[];
    Y=[];
    return
end

% Run over one mutiple of first guess for period to get better estimate.
total=30*Period;
Y0=X(end,:);

% Use events function to find period.
options=odeset('Events',@(t,X) PeriodFind(t,X,Y0),'RelTol',1e-8,'AbsTol',1e-10);
[~,X,te,~,~]=ode45(@(t,X) JR(t,X,P,A,B),[0 total],Y0,options);

% Retries using alternative events function if period is not found.
if length(te)<2
    disp('Period not found. Retrying...')
    
    % Events reports when y0 crosses from opposite direction.
    options=odeset('Events',@(t,X) PeriodFind2(t,X,Y0),'RelTol',1e-8,'AbsTol',1e-10);
    Y0=X(round((end/2)*(1+rand)),:);
    [~,X,te,~,~]=ode45(@(t,X) JR(t,X,P,A,B),[0 total],Y0,options);
end

% Calculates period.
Period=te(end)-te(end-1);

% Recomputes for mesh of uniform time points.
N=1000000;
dt=Period/N;
total=Period;
Y0=X(end,:);
options=odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,X]=ode45(@(t,X) JR(t,X,P,A,B),0:dt:total,Y0,options);

Y=X;

end

function JR=JR(t,X,P)

% Jansen-Rit ODEs

y=zeros(6,1);

y(1)=X(1);
y(2)=X(2);
y(3)=X(3);
y(4)=X(4);
y(5)=X(5);
y(6)=X(6);

JR=[y(4);y(5);y(6);
    P.A*P.a*sigm(y(2)-y(3))-2*P.a*y(4)-(P.a^2)*y(1);
    P.A*P.a*(P.P+P.C2*sigm(P.C1*y(1)))-2*P.a*y(5)-(P.a^2)*y(2);
    P.B*P.b*P.C4*sigm(P.C3*y(1))-2*P.b*y(6)-(P.b^2)*y(3)];
end

function sigm = sigm(v)

% Sigmoid

v0=6.0;
e0=2.5;
r=0.56;

sigm=(2*e0)./(1+exp(r*(v0-v)));

end

function [position,isterminal,direction] = PeriodFind(t,X,Y0)
position = X(1)-Y0(1); % Reports when y0 is the same as initial y0.
isterminal = 0;  % Halts integration when position=0
direction = 1;   % The zero is approached from positive direction
end

function [position,isterminal,direction] = PeriodFind(t,X,Y0)
position = X(1)-Y0(1); % Reports when y0 is the same as initial y0.
isterminal = 0;  % Halts integration when position=0
direction = -1;   % The zero is approached from negative direction
end