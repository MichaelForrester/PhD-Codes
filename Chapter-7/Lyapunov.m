function [R,LLE]=Lyapunov(S,f,P)

% Calculates LLE and rotation number for next-generation node with parameters P and sinusoidal input of
% magnitude S and frequency f.

% Set time increments.
T=200; dt=0.01; N=T/dt;

l=7; % Number of ODEs

% Random initial conditions
x0=[rand;rand;1/sqrt(2)*rand(2,1);rand;rand;0];

% Integrate forward 200 units
[~,x] = ode45(@(t,x) NextGen(t,x,P,S,freq),[0 1 200],x0);
% Begin at this point, hopefully near attractor!
x0=[x(end,1:end-1),0];

% Use events function to measure period.
options=odeset('Events',@(t,x) TheEventsFcn(t,x),'RelTol',1e-6,'AbsTol',1e-7);

% Integrate system
[t,x,te,ye,~] = ode45(@(t,x) NextGen(t,x,P,S,freq),0:dt:T,x0,options);

LLE=0; % Initialising.
w = [repmat(1/sqrt(6),6,1);0]; % Set random initial perturbation.
for k=1:N
    % calculate next point on trajectory.
    xn = x(k,:)';
    un = (eye(l)+NextGenJ(xn,P,S,freq)*dt)*w;
    % calculate stretching.
    a=norm(un);
    LLE=LLE+log(a);
    % renormalize.
    w=un/a;
end

LLE=LLE/T; % Exponent is given as average.

% Calculates rotation number if trajectory is stable
if LLE<0
    Periods=diff(te(abs(ye(:,1)-ye(end,1))<0.001))
    if size(Periods,1)>1 && Periods(end)-Periods(end-1)<0.01
        R=freq*Periods(end);
    end
else 
    R=0;
end
end

function [position,isterminal,direction] = TheEventsFcn(t,y)
    position = y(2); % The value that we want to be zero
    isterminal = 0;  % Don't halt integration 
    direction = 1;   % The zero is approached from positive direction
end