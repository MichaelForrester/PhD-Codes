function [F,Jac]=CBjac(S,P)
    
    Jac=zeros(14); % Preallocating Jacobian.
    F=zeros(14,1); % Preallocating function
    
    % Next-generation NMM variables
    xe=S(1); ye=S(2); xi=S(3); yi=S(4);
    gii=S(5); gie=S(6); gei=S(7); gee=S(8); gext=S(9);
    hii=S(10); hie=S(11); hei=S(12); hee=S(13); hext=S(14);
    
    % Input values into Jacobian
    Jac(1,1)=(1/P.taue)*(drfdx(xe,ye,P.ne,P.de)+drgdx(xe,ye,gei,P.vei)+drgdx(xe,ye,gee,P.vee)+drgdx(xe,ye,gext,P.vee));
    Jac(1,2)=(1/P.taue)*(drfdy(xe,ye,P.ne,P.de)+drgdy(xe,ye,gei,P.vei)+drgdy(xe,ye,gee,P.vee)+drgdy(xe,ye,gext,P.vee));
    Jac(2,1)=(1/P.taue)*(difdx(xe,ye,P.ne,P.de)+digdx(xe,ye,gei,P.vei)+digdx(xe,ye,gee,P.vee)+digdx(xe,ye,gext,P.vee));
    Jac(2,2)=(1/P.taue)*(difdy(xe,ye,P.ne,P.de)+digdy(xe,ye,gei,P.vei)+digdy(xe,ye,gee,P.vee)+digdy(xe,ye,gext,P.vee));

    Jac(3,3)=(1/P.taui)*(drfdx(xi,yi,P.ni,P.di)+drgdx(xi,yi,gii,P.vii)+drgdx(xi,yi,gie,P.vie));
    Jac(3,4)=(1/P.taui)*(drfdy(xi,yi,P.ni,P.di)+drgdy(xi,yi,gii,P.vii)+drgdy(xi,yi,gie,P.vie));
    Jac(4,3)=(1/P.taui)*(difdx(xi,yi,P.ni,P.di)+digdx(xi,yi,gii,P.vii)+digdx(xi,yi,gie,P.vie));
    Jac(4,4)=(1/P.taui)*(difdy(xi,yi,P.ni,P.di)+digdy(xi,yi,gii,P.vii)+digdy(xi,yi,gie,P.vie));

    Jac(1,7)=(1/P.taue)*drzdg(xe,ye,P.vei); Jac(1,8)=(1/P.taue)*drzdg(xe,ye,P.vee); Jac(1,9)=(1/P.taue)*drzdg(xe,ye,P.vee);
    Jac(2,7)=(1/P.taue)*dizdg(xe,ye,P.vei); Jac(2,8)=(1/P.taue)*dizdg(xe,ye,P.vee); Jac(2,9)=(1/P.taue)*dizdg(xe,ye,P.vee);
    Jac(3,5)=(1/P.taui)*drzdg(xi,yi,P.vii); Jac(3,6)=(1/P.taui)*drzdg(xi,yi,P.vie);
    Jac(4,5)=(1/P.taui)*dizdg(xi,yi,P.vii); Jac(4,6)=(1/P.taui)*dizdg(xi,yi,P.vie);

    Jac(5:9,10:14)=eye(5);

    Jac(10,5)=-P.aii^2; Jac(11,6)=-P.aie^2; Jac(12,7)=-P.aei^2;
    Jac(13,8)=-P.aee^2; Jac(14,9)=-P.aext^2;

    Jac(10,10)=-2*P.aii; Jac(11,11)=-2*P.aie; Jac(12,12)=-2*P.aei;
    Jac(13,13)=-2*P.aee; Jac(14,14)=-2*P.aext;

    Jac(10,3)=P.aii^2*P.kii*dfx(xi,yi,P.taui); Jac(10,4)=P.aii^2*P.kii*dfy(xi,yi,P.taui);
    Jac(11,1)=P.aie^2*P.kie*dfx(xe,ye,P.taue); Jac(11,2)=P.aie^2*P.kie*dfy(xe,ye,P.taue);
    Jac(12,3)=P.aei^2*P.kei*dfx(xi,yi,P.taui); Jac(12,4)=P.aei^2*P.kei*dfy(xi,yi,P.taui);
    Jac(13,1)=P.aee^2*P.kee*dfx(xe,ye,P.taue); Jac(13,2)=P.aee^2*P.kee*dfy(xe,ye,P.taue);
    Jac(14,1)=P.aext^2*P.e*dfx(xe,ye,P.taue); Jac(14,2)=P.aext^2*P.e*dfy(xe,ye,P.taue);
    
    % Compute ODE system RHS.
    F(1)=(1/P.taue)*(RF(xe,ye,P.ne,P.de)+RG(xe,ye,gei,P.vei)+RG(xe,ye,gee,P.vee)+RG(xe,ye,gext,P.vee));
    F(2)=(1/P.taue)*(IF(xe,ye,P.ne,P.de)+IG(xe,ye,gei,P.vei)+IG(xe,ye,gee,P.vee)+IG(xe,ye,gext,P.vee));
    F(3)=(1/P.taui)*(RF(xi,yi,P.ni,P.di)+RG(xi,yi,gii,P.vii)+RG(xi,yi,gie,P.vie));
    F(4)=(1/P.taui)*(IF(xi,yi,P.ni,P.di)+IG(xi,yi,gii,P.vii)+IG(xi,yi,gie,P.vie));

    F(5:9)=[hii;hie;hei;hee;hext];
    F(10)=P.aii^2*(P.kii*f(xi,yi,P.taui)-gii-2/P.aii*hii);
    F(11)=P.aie^2*(P.kie*f(xe,ye,P.taue)-gie-2/P.aie*hie);
    F(12)=P.aei^2*(P.kei*f(xi,yi,P.taui)-gei-2/P.aei*hei);
    F(13)=P.aee^2*(P.kee*f(xe,ye,P.taue)-gee-2/P.aee*hee);
    F(14)=P.aext^2*(P.e*f(xe,ye,P.taue)-gext-2/P.aext*hext); 
end

function X=RF(x,y,n,d)
    X=x*y-y-d*(x^2-y^2+2*x+1)/2+n*(-x*y-y);
end

function X=IF(x,y,n,d)
    X=(-x^2+y^2+2*x-1-d*(2*x*y+2*y)+n*(x^2-y^2+2*x+1))/2;
end

function X=RG(x,y,g,v)
    X=sum((-x*y-y)*v*g-(x^2-y^2-1)*g/2);
end

function X=IG(x,y,g,v)
    X=sum((x^2-y^2+2*x+1)*v*g/2-x*y*g);
end

function X=drfdx(x,y,n,d)
    X=y-d*(x+1)-y*n;
end

function X=drgdx(x,y,g,v)
    X=sum(-y*v*g-x*g);
end

function X=drfdy(x,y,n,d)
    X=x-1+d*y-n*(x+1);
end

function X=drgdy(x,y,g,v)
    X=sum(-(x+1)*v*g+y*g);
end

function X=difdx(x,y,n,d)
    X=-x+1-d*y+n*(x+1);
end

function X=digdx(x,y,g,v)
    X=sum((x+1)*v*g-y*g);
end

function X=difdy(x,y,n,d)
    X=y-d*(x+1)-n*y;
end

function X=digdy(x,y,g,v)
    X=sum(-y*v*g-x*g);
end

function X=drzdg(x,y,v)
    X=-v*y*(x+1)-(x^2-y^2-1)/2;
end

function X=dizdg(x,y,v)
    X=(x^2-y^2+2*x+1)*v/2-x*y;
end

function X = dfx(x,y,tau)
    X=-(1/(pi*tau))*2*(x.^2+2*x-y.^2+1)./((1+2*x+x.^2+y.^2).^2);
end

function X = dfy(x,y,tau)
    X=-(1/(pi*tau))*4*y.*(1+x)./((1+2*x+x.^2+y.^2).^2);
end

function f=f(x,y,tau)
    f=(1/(pi*tau))*(1-x.^2-y.^2)./(1+2*x+x.^2+y.^2);
end

