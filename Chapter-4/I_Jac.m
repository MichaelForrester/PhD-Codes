function J=I_Jac(pd,C,Y,u,Period)
    N=length(C); % Number of nodes.
    J=zeros(N,N); % Preallocating Jacobian
    
    if Period ~= 0
        
        % Calculates derivative of phase response curve.
        dt=Period/(length(Y(:,1))-1);
        udash=diff(u)/dt;
        udash(end+1)=(u(1)-u(end))/dt;
        for i=1:N
            % Calculates phase difference between node i and each connected node.
            phasedif=mod(pd(C(i,:)~=0)-pd(i),1);
            
            y=zeros(size(u,1),size(phasedif,1)); % Preallocating orbits.
            
            % Calculates phase-shifted orbits.
            for j=1:size(phasedif,1)
            y(:,j)=[Y(round(phasedif(j)*100000)+1:end);Y(1:round(phasedif(j)*100000))];
            end
            
            % Calculates derivative of phase interaction function.
            Idiff=-udash'*sigm(y)*dt/Period;
            
            % Calculates ith row of Jacobian at phase-locked state.
            J(i,C(i,:)~=0)=C(i,C(i,:)~=0).*Idiff;
            J(i,i)=-sum(J(i,:));
        end
    end
end

function sigm = sigm(v)

% Sigmoid

v0=6.0;
e0=2.5;
r=0.56;

sigm=(2*e0)./(1+exp(r*(v0-v)));

end