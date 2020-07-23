function I=PIF(pd,C,Y,u,Period)
    N=length(pd); % Number of nodes.
    I=zeros(N,1); % Preallocating
    L=length(Y(:,1)) % Number of time points
    dt=Period/(L-1); % Timeincrements
 
    % Computing phase interaction function
    if Period ~= 0
        for i=1:N
            phasedif=mod(pd(C(i,:)~=0)-pd(i),1);
            y=zeros(size(u,1),size(phasedif,1));
            for j=1:size(phasedif,1)
            y(:,j)=[Y(round(phasedif(j)*L)+1:end,2);Y(1:round(phasedif(j)*L),2)]-...
                [Y(round(phasedif(j)*L)+1:end,3);Y(1:round(phasedif(j)*L),3)];
            end
            I(i)=sum(C(i(C(i,:)~=0),:).*(u'*sigm(y)*dt/Period));
        end
    end
end

function sigm = sigm(v)

v0=6.0;
e0=2.5;
r=0.56;

sigm=(2*e0)./(1+exp(r*(v0-v)));

end