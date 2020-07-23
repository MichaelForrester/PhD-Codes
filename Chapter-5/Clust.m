function Clust=Clust(C,P)
    
    % Calculate multiplex clustering coefficient between SC and simulated FC.
    
    N=length(C); % Number of nodes.
    
    l=nnz(C); % Number of connections.
    
    tspan=[0,1000]; % Integration time limits.
    
    X0=rand(2*N,1); % Initial conditions for Wilson-Cowan network.
    
    [t,u]=WC_timeseties(P,C,X0,tspan); % Integrate system.
    
    % Continue if solution is oscillatory.
    [a,~] = findpeaks(u(:,1));
    [b,~] = findpeaks(-u(:,1));
    if max(a)-min(-b)>0.1
        
        % Compute FC
        U_trans=imag(hilbert(u));
        R=zeros(N);
        for n=1:N-1
            for m=n+1:N
                 R(n,m)=abs((1/length(t))*sum(exp(1i*(U_trans(:,m)-U_trans(:,n)))));
                 R(m,n)=R(n,m);
            end
        end
        
        % Compute clustering coefficient
        [~,S]=sort(R(:));
        Rthresh=zeros(N);
        Rthresh(S(end-l+1:end))=R(S(end-l+1:end));
        Rthresh(Rthresh~=0)=1; 
        num=C*(Rthresh.*(ones(N)-C))*C;
        den=C*(ones(N)-eye(N))*C-C^3;
        num=diag(num); den=diag(den);
        Clust=sum(num(den~=0)./den(den~=0))/N;
    else
        return;
    end
end