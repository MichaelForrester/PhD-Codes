# Produces time series for delayed system

# Loading packages
using DifferentialEquations
using MAT
using DSP
using SparseArrays
using LinearAlgebra

# Set parameters
const taui=5;  const taue=5;
const aee=0.6/taue; const aie=0.7/taui; const aei=1.2/taue; const aii=0.4/taui;
const kee=1.0*pi*taue; const kie=0.8*pi*taue; const kei=2*pi*taui; const kii=3*pi*taui;
const vee=6; const vie=10; const vei=-10; const vii=-10;
const de = 0.5; const di = 0.5;
const kc=10.0;
const ni=-40; const ne=20;
cv=0.01; od=0;

# Load SC matrix
file = matopen("C.mat");
C=0.1*read(file, "C");
close(file);

const Nn=size(C,1); # Number of nodes

# Adding self-connections
const C=C+kee*Diagonal((ones(Nn,Nn)))

# Loading initial conditions
file = matopen("h.mat");
const IC=read(file, "h");
close(file);

const Nc=sum(!iszero,C); # Storing node degrees# Number of connections

# Storing node degrees
Cn_temp=zeros(size(C));
Cn_temp[findall(!iszero,C)]=ones(size(findall(!iszero,C)));
const Cn=append!([0],sum(Cn_temp,dims=2));

# Storing number of ODEs for each node
Vn=zeros(1,Nn+1);
const Vn[2:end]=cumsum(2*transpose(Cn[2:end]).+10,dims=2);
const C_lags=cumsum(Cn-append!([0],repeat([1],Nn,1)),dims=1);

# History function
h(p, t; idxs = nothing)=1;

# Defining functions
f(x,y,tau)=(1/(pi*tau))*(1 .-x.^2-y.^2)./(1 .+2*x+x.^2+y.^2)
RF(x,y,n,d)=x.*y-y-d*(x.^2-y.^2+2*x.+1)/2+n*(-x.*y-y);
IF(x,y,n,d)=(-x.^2+y.^2+2*x.-1-d*(2*x.*y+2*y)+n*(x.^2-y.^2+2*x.+1))/2;
RG(x,y,g,v)=sum((-x.*y-y)*v.*g-(x.^2-y.^2 .-1).*g/2);
IG(x,y,g,v)=sum((x.^2-y.^2+2*x.+1)*v.*g/2-x.*y.*g);

# Defining next-generation neural mass model
function CB_TMS(dydt,S,h,p,t)
    for n=1:Nn
        o=trunc(Int,Vn[n]);
        co=trunc(Int,Cn[n+1]);
        Fn=zeros(trunc(Int, Vn[n+1]-Vn[n]),1);
        xe=[S[o+1]]; ye=[S[o+2]]; xi=[S[o+3]]; yi=[S[o+4]];
        gii=[S[o+5]]; gie=[S[o+6]]; gei=[S[o+7]]; gee=S[o+8:o+7+co];
        hii=[S[o+8+co]]; hie=[S[o+9+co]]; hei=[S[o+10+co]]; hee=S[o+11+co:o+10+2*co];

        Fn[1]=(1/taue)*(RF(xe,ye,ne,de)[1]+RG(xe,ye,gei,vei)+RG(xe,ye,gee,vee));
        Fn[2]=(1/taue)*(IF(xe,ye,ne,de)[1]+IG(xe,ye,gei,vei)+IG(xe,ye,gee,vee));
        Fn[3]=(1/taui)*(RF(xi,yi,ni,di)[1]+RG(xi,yi,gii,vii)+RG(xi,yi,gie,vie));
        Fn[4]=(1/taui)*(IF(xi,yi,ni,di)[1]+IG(xi,yi,gii,vii)+IG(xi,yi,gie,vie));

        Fn[5:7+co]=[hii;hie;hei;hee];
        Fn[8+co]=(aii^2*(kii*f(xi,yi,taui)-gii-2/aii*hii))[1];
        Fn[9+co]=(aie^2*(kie*f(xe,ye,taue)-gie-2/aie*hie))[1];
        Fn[10+co]=(aei^2*(kei*f(xi,yi,taui)-gei-2/aei*hei))[1];
        f_del=zeros(co,1);
        f_del[findall(C[n,C[n,:].!=0].==kee)].=f(xe,ye,taue);
        del_num=findall(iszero,f_del);
        for m=1:co-1
            f_del[del_num[m]]=f(h(p,t-lags[C_lags[n]+m];idxs=round.(Int,Vn[findall(!iszero,C[n,:])[del_num[m]]])+1),h(p,t-lags[C_lags[n]+m];idxs=round.(Int,Vn[findall(!iszero,C[n,:])[del_num[m]]])+2),taue)[1];
        end

        Fn[11+co:10+2*co]=aee^2*(C[findall(!iszero,C[:,n]),n].*f_del-gee-2/aee*hee);

        dydt[o+1:trunc(Int,Vn[n+1])]=Fn;
    end
end

# Set delays
file = matopen("lags.mat");
lag=round.(read(file, "lags")/cv);
close(file);
const lags= lag+repeat([od],length(lag),1);
lag=nothing;

# Integrate system
tspan = (0.0,600.0);
u0 = IC;
prob = DDEProblem(CB_TMS,u0,h,tspan; constant_lags=lags);
alg = MethodOfSteps(Tsit5());
sol=solve(prob, maxiters = 1e8, alg, progress=true, abstol = 1e-8, reltol = 1e-8, save_idxs = append!(round.(Int,Vn[1:end-1].+1)[:,1],round.(Int,Vn[1:end-1].+2)[:,1]));

# Save time series for firing rate of excitatory populations
U=zeros(size(sol.u)[1],2*Nn);
for n=1:size(sol.u)[1]
    U[n,:]=sol.u[n];
end
U=f(U[:,1:Nn],U[:,Nn+1:end],taue);

# Compute FC
U_trans=imag(hilbert(U));
R=zeros(Nn,Nn);
for n=1:Nn-1
    for m=n+1:Nn
        R[n,m]=abs((1/size(U_trans)[1]*sum(exp.(im*(U_trans[:,m]-U_trans[:,n])))));
        R[m,n]=R[n,m];
    end
end

# Save time series and FC
file = matopen("timeseries.mat", "w");
write(file, "U", U);
close(file);

file = matopen("timeseries_t.mat", "w");
write(file, "t", sol.t);
close(file);

file = matopen("FC.mat", "w");
write(file, "R", R);
close(file);
