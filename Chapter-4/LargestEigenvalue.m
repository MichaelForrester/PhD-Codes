function LE=LargestEigenvalue(C,P)

N=length(C,1); % Number of nodes.
T = 100; % Total time for time-series
U0=rand(1,6); % Initial conditions for phase response curve evaluation.
yi=rand(1,6); % Initial conditions for Jansen-Rit time-series

% Compute phase response curve U, orbit Y and period of oscillation.
[Period,~,Y,U] = PRC(P,U0,yi);
u=U(:,5); %  Only the 5th variable is dependent on interactions.

% Loop ensures we have N phase differences.
ie=0;
while length(ie)<N
    % Set initial conditions chosen randomly from orbit.
    y_i=Y(round(size(Y,1)*rand(1,N)),:); y_i=y_i(:);
    
    % Run until total time (for convergence to phase-locked state).
    [~,y]=ode45(@(t,y) JansenRit(t,y,P,C,N), [0,T/2,T], y_i);
    
    % Reset initial conditions.
    y_i=y(end,:);
    
    % Use function PDFind to find define phases (relative to 1st node).
    opts = odeset('Events',@(t,X) PDFind(t,X,y_i(1:N)),'RelTol',1e-6,'AbsTol',1e-8);
    
    % Reports times, te, that nodes, ie, correspond with the initial phase of 1st node.
    [~,~,te,~,ie] =ode45(@(t,y) JansenRit(t,y,P,C,N), [0 5*Period], y_i,opts);
end

% Defines variable determining interaction.
H=Y(:,2)-Y(:,3);

% Calculates set of phase differences.
PD=PDForm(ie,te,N,Period);

% Calculates eigenvalues of Jacobian at phase-locked state.
vals=eig(I_Jac(PD,P,H,u,Period));

% Reports eigenvalue with maximum real part.
LE=max(real(vals));

end

function [position,isterminal,direction] = PDFind(t,X,Y0)
    % Evaluates when the variable y0 each node crosses a reference point in the orbit,
    % defined to be the initial value for the first node. 

    % Calculates difference between variable y0 for each node an initial y0 of 1st node. 
    position = X(1:78)-Y0(1);
    
    % Halts integration when condition is met for any node.
    isterminal = zeros(78,1);
    
    % The zero must be approached from positive direction.
    direction = ones(78,1);
end

function PD=PDForm(ie,te,N,Period)
    % Calculates set of phase differences.

    % Takes final crossing times and node labels for each node.
    t=te(end-N+1:end); [~,I]=sort(ie(end-N+1:end));
    
    PD=t(I); % Sorts phase differences from 1 to N.
    PD=PD-PD(1); % Shifts phase differences relative to 1st node.
    PD=mod(PD,Period); % Ensures all phase differences are within [0,Period].
    PD=PD/Period; % Normalising.
end
