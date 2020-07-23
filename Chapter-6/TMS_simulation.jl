# Simulation of TMS in human cortical network

# Loading packages
using DifferentialEquations
using MAT
using DSP
using SparseArrays
using LinearAlgebra

# Setting parameters
const aee=1; const aie=1.4; const aei=0.7; const aii=0.4;
const kee=1.5; const kie=1; const kei=2; const kii=3;
const vee=10; const vie=8; const vei=-8; const vii=-12;
const de = 0.5; const di = 0.5;
const ni0=-47; const ne0=20;

const target=1; # Target region

# Load SC matrix
file = matopen("C.mat");
C=0.1*read(file, "C");
close(file);

const Nn=size(C,1); # Number of nodes

# Adding self-connections
const C=C+kee*Diagonal((ones(Nn,Nn)))

# Load initial conditions
file = matopen("init.mat");
const IC=read(file, "S");
close(file);

const Nc=sum(!iszero,C); # Number of connections

# Storing node degrees
Cn_temp=zeros(size(C));
Cn_temp[findall(!iszero,C)]=ones(size(findall(!iszero,C)));
const Cn=sum(Cn_temp,dims=1);

# Storing number of ODEs for each node
Vn=zeros(1,Nn+1);
const Vn[2:end]=cumsum(2*Cn+repeat([10],1,Nn),dims=2);

# Defining functions
f(x,y)=(repeat([1],length(x),1)-x.^2-y.^2)./(repeat([1],length(x),1)+2*x+x.^2+y.^2)
RF(x,y,n,d)=x.*y-y-d*(x.^2-y.^2+2*x+repeat([1],length(x),1))/2+n*(-x.*y-y);
IF(x,y,n,d)=(-x.^2+y.^2+2*x-repeat([1],length(x),1)-d*(2*x.*y+2*y)+n*(x.^2-y.^2+2*x+repeat([1],length(x),1)))/2;
RG(x,y,g,v)=sum((-x.*y-y)*v.*g-(x.^2-y.^2-repeat([1],length(x),1)).*g/2);
IG(x,y,g,v)=sum((x.^2-y.^2+2*x+repeat([1],length(x),1))*v.*g/2-x.*y.*g);

# TMS pulse function
function TMS_pulse(t,freq,peak)
    w=20000; tau=0.00008;
    ts = mod(t,1/freq);
    if t>100 && t<150
        pulse=peak*sin(w*ts)*exp(-ts/tau);
        return pulse;
    else
        return 0;
    end
end

# Next-generation neural mass model
function CB_TMS(dydt,S,p,t)

    for n=1:Nn

        if n==target
            ne=ne0+TMS_pulse(t,20.0,100.0);
            ni=ni0+TMS_pulse(t,20.0,100.0);
        else
            ne=ne0;
            ni=ni0;
        end

        o=trunc(Int,Vn[n]);
        co=trunc(Int,Cn[n]);
        Fn=zeros(trunc(Int, Vn[n+1]-Vn[n]),1);
        xe=[S[o+1]]; ye=[S[o+2]]; xi=[S[o+3]]; yi=[S[o+4]];
        gii=[S[o+5]]; gie=[S[o+6]]; gei=[S[o+7]]; gee=S[o+8:o+7+co];
        hii=[S[o+8+co]]; hie=[S[o+9+co]]; hei=[S[o+10+co]]; hee=S[o+11+co:o+10+2*co];

        Fn[1]=RF(xe,ye,ne,de)[1]+RG(xe,ye,gei,vei)+RG(xe,ye,gee,vee);
        Fn[2]=IF(xe,ye,ne,de)[1]+IG(xe,ye,gei,vei)+IG(xe,ye,gee,vee);
        Fn[3]=RF(xi,yi,ni,di)[1]+RG(xi,yi,gii,vii)+RG(xi,yi,gie,vie);
        Fn[4]=IF(xi,yi,ni,di)[1]+IG(xi,yi,gii,vii)+IG(xi,yi,gie,vie);

        Fn[5:7+co]=[hii;hie;hei;hee];
        Fn[8+co]=dropdims(aii^2*(kii*f(xi,yi)-gii-2/aii*hii),dims=1)[1];
        Fn[9+co]=dropdims(aie^2*(kie*f(xe,ye)-gie-2/aie*hie),dims=1)[1];
        Fn[10+co]=dropdims(aei^2*(kei*f(xi,yi)-gei-2/aei*hei),dims=1)[1];
        Fn[11+co:10+2*co]=aee^2*(C[findall(!iszero,C[:,n]),n].*f(S[round.(Int,Vn[findall(!iszero,C[n,:])])+repeat([1],co,1)],S[round.(Int,Vn[findall(!iszero,C[n,:])])+repeat([2],co,1)])-gee-2/aee*hee);

        dydt[o+1:trunc(Int,Vn[n+1])]=Fn;
    end
end

# Integrating system
tspan = (0.0,300.0);
u0 = IC;
prob = ODEProblem(CB_TMS,u0,tspan);
alg = MethodOfSteps(Tsit5());
sol=solve(prob, maxiters = 1e8, alg, progress=true);

# Save solution and time.
file = matopen("steadystate.mat", "w");
write(file, "sol", sol.u);
close(file);

file = matopen("steadystate_time.mat", "w");
write(file, "t", sol.t);
close(file);
