function f=NextGen(t,X,P,S,freq)

% Next-generation NMM RHS

g=X(1); h=X(2); x=X(3); y=X(4); A=X(5); dA=X(6); tau=X(7);

f=zeros(7,1);

f(1)=h;
f(2)=P.a^2*(P.k*fr(x,y)-g-2/P.a*h);
f(3)=RF(x,y,P.n+A,P.d)+RG(x,y,g,P.v);
f(4)=IF(x,y,P.n+A,P.d)+IG(x,y,g,P.v);
f(5)=dA;
f(6)=P.ad^2*(I(tau,S,freq)-A-2/P.ad*dA);
f(7)=1;

end

    function X=RF(x,y,n,d)
        X=x*y-y-d*(x^2-y^2+2*x+1)/2+n*(-x*y-y);
    end

    function X=IF(x,y,n,d)
        X=(-x^2+y^2+2*x-1-d*(2*x*y+2*y)+n*(x^2-y^2+2*x+1))/2;
    end

    function X=RG(x,y,g,v)
        X=(-x*y-y)*v*g-(x^2-y^2-1)*g/2;
    end

    function X=IG(x,y,g,v)
        X=(x^2-y^2+2*x+1)*v*g/2-x*y*g;
    end

    function f=fr(x,y)
        f=(1/pi)*(1-x.^2-y.^2)./(1+2*x+x.^2+y.^2);
    end

    function X=I(t,S,f)
        X=S/2*(1+sin(f*2*pi*t));
    end

%     function X=I(t,S,f)
%         X=S.*(mod(t,f)<0.1);
%     end