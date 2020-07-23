function [nodel_spectra,del_spectra]=CBdelplot(C,D,y0,P,s)
    % Calculates the spectra for an undelayed network of next-generation neural masses with paramters set in P, as well as spectra for the same network with conduction delays with distance matrix D and conduction speed s.
    
    N=length(C); % Number of nodes.
    M=length(y0); % Number of ODEs for each node.
    
    % Calculating delay matrix (distance/speed).
    tau=D/s;
    
    C0=P.e*C; % Scaling coupling.
    
    [EVec,EVal]=eig(C0); % Eigenvectors and eigenvalues of SC matrix.
    
    % Defining solver options.
    options = optimoptions(@fsolve,'Display','iter',...
        'Algorithm','trust-region-dogleg',...
        'SpecifyObjectiveGradient',true,'PrecondBandWidth',0);
    
    % Solving neural-mass network using next-generation model and its jacobian defined in NextGen, with initial conditions y0.
    [~,~,~,~,jacobian] = fsolve(@(y) NextGen(y,P),y0,options);
    
    jacobian=full(jacobian); % Conversion from sparse to full Jacobian.
    
    % Setting DF (from equation 7.8).
    DF=jacobian; DF(M,1:2)=0;
    
    % Setting DG (from equation 7.8).
    DG=zeros(M);
    DG(M,1)=P.aext^2*dfx(S(1),S(2),P.taue);
    DG(M,2)=P.aext^2*dfy(S(1),S(2),P.taue);
    
    % Redefining options
    options = optimoptions(@fsolve,'Display','iter',...
        'Algorithm','trust-region-dogleg',...
        'SpecifyObjectiveGradient',false,'PrecondBandWidth',0);
    
    nodel_spectra=zeros(M*N,1); % Initialising spectra for undelayed network.
    del_spectra=[]; % Initialising spectra for delayed network.
    
    % Loop over each eigenmode n.
    for n=1:N
        
        % Spectra for nth eigenmode of undelayed system.
        nodel_spectra((n-1)*M+1:n*M)=eig(DF+EVal(n)*DG);
        
        % Set upper and lower limits for real and imaginary eigenvalues and make arrays with 1000 uniform intervals.
        u_lim=[u_min u_max]; v_lim=[v_min v_max];
        u=linspace(u_lim(1),u_lim(2),1001);
        v=linspace(v_lim(1),v_lim(2),1001);
        
        % Initialise plane in selected domain of real/imaginary values.
        UVplane=zeros(length(v),length(u));
        
        % Calculate absolute value of characteristic determinant for each complex value.
        for U=1:length(u)
            for V=1:length(v)
                ig=u(U)+1i*v(V);
                UVplane(V,U)=abs(del_det(l,C0,tau,DF,DG,EVec(:,n)));
            end
        end
        
        % Find local minima in UVplane
        [X,Y]=meshgrid(u,v);
        ix = find(imregionalmin(UVplane));
        
        % Setting initial conditions for solver
        X0=X(ix)+1i*Y(ix);
        
        % Solving for each eigenvalue
        for m=1:length(X0)
            x0=X0(m);
            [l,~,ef,~] = fsolve(@(l) del_det(l,C0,tau,DF,DG,EVec(:,n)),x0,options);
            if ef~=-2 && ef~=-3
                del_spectra=[del_spectra;real(l),imag(l),n];
            end
        end
    end
end

% Defining df/dx and df/dy, where f is firing rate of excitatory population, x and y are the real and imaginary Kuramoto order paramters, respectively.
function X = dfx(x,y,tau)
    X=-(1/(pi*tau))*2*(x.^2+2*x-y.^2+1)./((1+2*x+x.^2+y.^2).^2);
end

function X = dfy(x,y,tau)
    X=-(1/(pi*tau))*4*y.*(1+x)./((1+2*x+x.^2+y.^2).^2);
end

% Characteristic determinant.
function X=del_det(l,C,tau,DF,DG,V)
    X=det(eye(length(DF))*l-(DF+sum(sum(C.*exp(-l*tau).*(V*V')))*DG));
end