% ODEsolve Weakly-Coupled

function dydt=WeakCoup(t,y,P,Y,u,Period,F)
    pd=mod((y-y(1))/Period,1); % Normalise phase differences
    dydt=F+P.A*P.a*P.ep*I(pd,C,Y,u,Period); % RHS of oscillator network
end