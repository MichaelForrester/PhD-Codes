function Jaccard=JansenRit(C,P)

%Function to calculate the Jaccard coefficient for the SC matrix C and the FC matrix
%resulting from the MPC between each node's timeseries (transformed by the Hilbert function).
%Timeseries dynamics are governed by the Jansen-Rit model with parameters P.

% Inputs: C - Connectivity matrix where C(i,j) is coupling strength from node j -> i.
%         P - Structure containing Jansen-Rit paramteters.
% Output: Jaccard - Jaccard coefficient between SC and simulated FC.
.
N=length(C); % Number of nodes.

% Setting time increments.
T = 500; dt = 0.001; N = T/dt;

% Preallocating variables.
y0=rand(N,1);
y1=rand(N,1);
y2=rand(N,1);
y3=rand(N,1);
y4=rand(N,1);
y5=rand(N,1);

% Preallocating saved variables.
y0save = zeros(N,1);
y1save = zeros(N,1);
y2save = zeros(N,1);
y3save = zeros(N,1);
y4save = zeros(N,1);
y5save = zeros(N,1);
ysave  = zeros(N,N+1);

% Preallocating space for timeseries of variable y.
ysave(:,1) = y1-y2;

% Defining noise.
dW = normrnd(0,0.1,F,N-1);

% Implementing method.
for j = 1:N
        y0save = y0 + y3*dt;
        y1save = y1 + y4*dt;
        y2save = y2 + y5*dt;
        y3save = y3 + ...
            (P.A*P.a*sigm(y1-y2)-2*P.a*y3-(P.a^2)*y0)*dt;
        y4save = y4 + ...
            (P.A*P.a*(P.p.P+P.e*(C*f(y1-y2))+dW(:,N-1)+P.C2*f(P.C1*y0))-2*P.a*y4-(P.a^2)*y1)*dt;
        y5save = y5 + ...
            (P.B*P.b*P.C4*f(P.C3*y0)-2*P.b*y5-(P.b^2)*y2)*dt;
        
        y0 = y0save; y1 = y1save; y2 = y2save;
        y3 = y3save; y4 = y4save; y5 = y5save;
        ysave(:,j+1)=y1-y2;
end

ysave(:,1:100000)=[]; % Removing inital transients.
U_trans = angle(hilbert(ysave')); % Implementing Hilbert tranform.
FC = zeros(N); % Preallocating FC.

% Calculating FC matrix.
for f = 1:F
    FC(f,f+1:F) = abs((1/size(U_trans,1))*sum(exp(1i*(unwrap(U_trans(:,f+1:F))-repmat(unwrap(U_trans(:,f)),1,F-f)))));
    FC(f+1:F,f)=FC(f,f+1:F);
end

% Thresholding and binarising FC matrix.
R_bin=zeros(F);
[~,I]=sort(FC(:));
R_bin(I(end-cons+1:end))=1;

% Calculating Jaccard coefficient.
Jaccard = sum(sum(Connect_bin.*R_bin))/sum(sum(Connect_bin==1|R_bin==1));
end
