function [FCjr,FCwc]=FC(C,P)

% Compute FC matrices for Jansen-Rit network and weakly-coupled reduction.

N=length(C); % Number of nodes.

% Set time increments.
T = 1000; dt=0.001; tspan=0:dt:T;

% Compute phase response function
yi=rand(1,6);
U0=rand(1,6);
[Period,~,Y,U] = PRC(P,U0,yi);
u=U(:,5);

Tn=length(Y(:,1)); % Number of time increments of Y.

freq=1; % Normalised frequency

L=size(tspan,2); 

FCwc=zeros(N);
FCjr=zeros(N);
        
    p_i=Period*rand(N,1); % Initial conditions for phase.

    y_i=zeros(N*6,1); % Initial conditions for Jansen-Rit system.
    for n=1:F
        y_i((0:5)*N+n)=Y(round((p_i(n)/Period)*Tn),:);
    end

    % Integrate weakly-coupled oscillator system.
    [~,y] = ode45(@(t,y) WeakCoup(t,y,A,B,Y,u,Period,freq,ep), tspan, p_i);
    phase=mod(y,Period);

    % Compute MPC for weakly-coupled oscillator time series.
    U_trans=phase*2*pi/Period;
    R = zeros(F);
    for f = 1:F
        R(f,f+1:F) = abs((1/L)*sum(exp(1i*(U_trans(:,f+1:F)-repmat(U_trans(:,f),1,F-f)))));
        R(f+1:F,f)=R(f,f+1:F);
    end
    FCwc=R;
    
    % Integrate  Jasen-Rit system
    [~,y] = ode45(@(t,y) JansenRit(t,y,A,B,P,ep,Connect,F), tspan, y_i);
    
    % Compute MPC for Jansen-Rit time series.
    y=y(:,F+1:2*F)-y(:,2*F+1:3*F);
    y=y-mean(y);
    U_trans=angle(hilbert(y));
    for f = 1:F
        R(f,f+1:F) = abs((1/L)*sum(exp(1i*(unwrap(U_trans(:,f+1:F))-repmat(unwrap(U_trans(:,f)),1,F-f)))));
        R(f+1:F,f)=R(f,f+1:F);
    end
    FCjr=R;

end